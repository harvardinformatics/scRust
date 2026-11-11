//! JSON in the layout of the reference scripts' `to_json` (`validation/R/lib/io.R`), so that a
//! Rust `summary.json` reads like R's and its `result` section can be compared text to text:
//! objects one member per line with two-space indentation, vectors of scalars on one line
//! (`[1, 2]`), lists of objects or arrays one item per line, doubles as `%.17g` (non-finite as
//! `null`). A small parser reads what R writes back, keeping each number's text.

use crate::text::fmt_g;
use std::fmt::Write as _;

/// A JSON value.
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    /// An integer (R's integer type), written as is.
    Int(i64),
    /// A double, written as `%.17g`; NA, NaN and infinities as `null`.
    Num(f64),
    Str(String),
    /// A vector of scalars, written on one line.
    Vector(Vec<Json>),
    /// An unnamed list, one item per line.
    List(Vec<Json>),
    /// A named list, one member per line, in order.
    Object(Vec<(String, Json)>),
}

impl Json {
    pub fn obj(members: Vec<(&str, Json)>) -> Json {
        Json::Object(members.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
    }

    pub fn strs<S: AsRef<str>>(v: &[S]) -> Json {
        Json::Vector(v.iter().map(|s| Json::Str(s.as_ref().to_string())).collect())
    }

    pub fn nums(v: &[f64]) -> Json {
        Json::Vector(v.iter().map(|&x| Json::Num(x)).collect())
    }

    /// The member `key` of an object.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(m) => m.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(x) => Some(*x),
            Json::Int(i) => Some(*i as f64),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    /// The text R's `to_json` writes for this value at `indent`.
    pub fn to_text(&self, indent: usize) -> String {
        let pad = "  ".repeat(indent);
        let pad1 = "  ".repeat(indent + 1);
        match self {
            Json::Null => "null".into(),
            Json::Bool(b) => b.to_string(),
            Json::Int(i) => i.to_string(),
            Json::Num(x) => {
                if x.is_finite() {
                    fmt_g(*x, 17)
                } else {
                    "null".into()
                }
            }
            Json::Str(s) => escape(s),
            Json::Vector(v) => {
                let items: Vec<String> = v.iter().map(|x| x.to_text(indent + 1)).collect();
                format!("[{}]", items.join(", "))
            }
            Json::List(v) => {
                if v.is_empty() {
                    return "[]".into();
                }
                let items: Vec<String> = v.iter().map(|x| format!("{pad1}{}", x.to_text(indent + 1))).collect();
                format!("[\n{}\n{pad}]", items.join(",\n"))
            }
            Json::Object(m) => {
                if m.is_empty() {
                    return "{}".into();
                }
                let items: Vec<String> =
                    m.iter().map(|(k, v)| format!("{pad1}{}: {}", escape(k), v.to_text(indent + 1))).collect();
                format!("{{\n{}\n{pad}}}", items.join(",\n"))
            }
        }
    }

    /// Parses JSON text. Numbers with a fraction or exponent become [`Json::Num`] (Rust's
    /// parser is correctly rounded, so `%.17g` text gives back the exact double); others
    /// [`Json::Int`]. Arrays of scalars become [`Json::Vector`], other arrays [`Json::List`].
    pub fn parse(text: &str) -> Result<Json, String> {
        let mut p = Parser { s: text.as_bytes(), i: 0 };
        let v = p.value()?;
        p.ws();
        if p.i != p.s.len() {
            return Err(format!("trailing text at byte {}", p.i));
        }
        Ok(v)
    }
}

fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }

    fn expect(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.s.get(self.i) == Some(&c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected '{}' at byte {}", c as char, self.i))
        }
    }

    fn value(&mut self) -> Result<Json, String> {
        self.ws();
        match self.s.get(self.i) {
            Some(b'{') => {
                self.i += 1;
                let mut m = Vec::new();
                self.ws();
                if self.s.get(self.i) == Some(&b'}') {
                    self.i += 1;
                    return Ok(Json::Object(m));
                }
                loop {
                    self.ws();
                    let k = self.string()?;
                    self.expect(b':')?;
                    m.push((k, self.value()?));
                    self.ws();
                    match self.s.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b'}') => {
                            self.i += 1;
                            return Ok(Json::Object(m));
                        }
                        _ => return Err(format!("expected ',' or '}}' at byte {}", self.i)),
                    }
                }
            }
            Some(b'[') => {
                self.i += 1;
                let mut v = Vec::new();
                self.ws();
                if self.s.get(self.i) == Some(&b']') {
                    self.i += 1;
                    return Ok(Json::Vector(v));
                }
                loop {
                    v.push(self.value()?);
                    self.ws();
                    match self.s.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b']') => {
                            self.i += 1;
                            let scalars = v.iter().all(|x| !matches!(x, Json::Vector(_) | Json::List(_) | Json::Object(_)));
                            return Ok(if scalars { Json::Vector(v) } else { Json::List(v) });
                        }
                        _ => return Err(format!("expected ',' or ']' at byte {}", self.i)),
                    }
                }
            }
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b't') if self.s[self.i..].starts_with(b"true") => {
                self.i += 4;
                Ok(Json::Bool(true))
            }
            Some(b'f') if self.s[self.i..].starts_with(b"false") => {
                self.i += 5;
                Ok(Json::Bool(false))
            }
            Some(b'n') if self.s[self.i..].starts_with(b"null") => {
                self.i += 4;
                Ok(Json::Null)
            }
            Some(_) => {
                let start = self.i;
                while self.i < self.s.len() && matches!(self.s[self.i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.s[start..self.i]).unwrap();
                if t.is_empty() {
                    return Err(format!("unexpected character at byte {start}"));
                }
                if t.contains(['.', 'e', 'E']) {
                    t.parse().map(Json::Num).map_err(|_| format!("bad number {t}"))
                } else {
                    t.parse().map(Json::Int).map_err(|_| format!("bad number {t}"))
                }
            }
            None => Err("unexpected end of JSON".into()),
        }
    }

    fn string(&mut self) -> Result<String, String> {
        if self.s.get(self.i) != Some(&b'"') {
            return Err(format!("expected a string at byte {}", self.i));
        }
        self.i += 1;
        let mut out = String::new();
        loop {
            let start = self.i;
            while self.i < self.s.len() && self.s[self.i] != b'"' && self.s[self.i] != b'\\' {
                self.i += 1;
            }
            out.push_str(std::str::from_utf8(&self.s[start..self.i]).map_err(|_| "invalid UTF-8")?);
            match self.s.get(self.i) {
                Some(b'"') => {
                    self.i += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    let c = *self.s.get(self.i + 1).ok_or("bad escape")?;
                    self.i += 2;
                    match c {
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'u' => {
                            let h = std::str::from_utf8(self.s.get(self.i..self.i + 4).ok_or("bad \\u")?).unwrap();
                            let cp = u32::from_str_radix(h, 16).map_err(|_| "bad \\u")?;
                            out.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                            self.i += 4;
                        }
                        c => out.push(c as char),
                    }
                }
                _ => return Err("unterminated string".into()),
            }
        }
    }
}

/// Writes `value` followed by a newline, as R's `writeLines(to_json(x))`.
pub fn to_file_text(value: &Json) -> String {
    let mut s = value.to_text(0);
    let _ = writeln!(s);
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_layout_round_trip() {
        let v = Json::obj(vec![
            ("result", Json::obj(vec![
                ("n_cells", Json::Int(300)),
                ("pK", Json::Num(0.29)),
                ("labels", Json::strs(&["0.04", "0.05"])),
                ("k", Json::Num(116.0)),
                ("none", Json::Null),
            ])),
            ("flag", Json::Bool(true)),
        ]);
        let text = to_file_text(&v);
        assert_eq!(
            text,
            "{\n  \"result\": {\n    \"n_cells\": 300,\n    \"pK\": 0.28999999999999998,\n    \
             \"labels\": [\"0.04\", \"0.05\"],\n    \"k\": 116,\n    \"none\": null\n  },\n  \"flag\": true\n}\n"
        );
        let back = Json::parse(&text).unwrap();
        assert_eq!(back.get("result").unwrap().get("pK").unwrap().as_f64(), Some(0.29));
        assert_eq!(back.get("result").unwrap().get("k"), Some(&Json::Int(116)));
    }
}

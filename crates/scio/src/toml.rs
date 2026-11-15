//! A small TOML reader (TOML 1.0, standard library only), for parameter files
//! (`scrust cluster --params`, docs/PLAN.md C3) and the repository's own TOML files
//! (`validation/tolerances.toml`, `validation/manifest.toml`).
//!
//! It reads tables (`[a.b]`), arrays of tables (`[[a]]`), dotted keys, bare and quoted keys,
//! basic, literal and multi-line strings, integers (decimal, hex, octal, binary, with `_`),
//! floats (with exponents, `inf`, `nan`), booleans, arrays (over several lines, with comments
//! and a trailing comma) and inline tables. Dates and times are refused. A key defined twice, or a
//! table defined twice, is an error, as TOML requires.
//!
//! **Numbers keep their text.** A float holds the value Rust's correctly rounded parser gives and
//! the text as written, so a caller that must read a number as R does (`as.numeric`, R's own
//! decimal parser, which is not always correctly rounded) can parse the text again.

use std::fmt;

/// A TOML value. Tables keep their keys in the order of the file.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    String(String),
    /// An integer, and its text as written (without `_`).
    Integer(i64, String),
    /// A float, and its text as written (without `_`).
    Float(f64, String),
    Boolean(bool),
    Array(Vec<Value>),
    Table(Table),
}

/// A table: keys in file order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Table {
    pub entries: Vec<(String, Value)>,
    /// Defined by a header (`[a]`) or inline (`{...}`), as opposed to implied by a dotted key or
    /// a sub-table's header; TOML forbids defining a table twice.
    defined: bool,
    /// An inline table or an array's table: closed to later additions.
    sealed: bool,
}

impl Table {
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.entries.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    fn get_mut(&mut self, key: &str) -> Option<&mut Value> {
        self.entries.iter_mut().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// The value at a dotted path, `a.b.c`.
    pub fn path(&self, path: &str) -> Option<&Value> {
        let mut parts = path.split('.');
        let mut v = self.get(parts.next()?)?;
        for p in parts {
            v = match v {
                Value::Table(t) => t.get(p)?,
                _ => return None,
            };
        }
        Some(v)
    }

    pub fn keys(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|(k, _)| k.as_str())
    }
}

impl Value {
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::String(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Integer(i, _) => Some(*i),
            _ => None,
        }
    }

    /// A number as a double (integers converted).
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Integer(i, _) => Some(*i as f64),
            Value::Float(x, _) => Some(*x),
            _ => None,
        }
    }

    /// A number's text as written (decimal integers and floats), for callers that parse it again.
    pub fn number_text(&self) -> Option<&str> {
        match self {
            Value::Integer(_, t) | Value::Float(_, t) => Some(t),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Boolean(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[Value]> {
        match self {
            Value::Array(a) => Some(a),
            _ => None,
        }
    }

    pub fn as_table(&self) -> Option<&Table> {
        match self {
            Value::Table(t) => Some(t),
            _ => None,
        }
    }

    /// What kind of value this is, for messages.
    pub fn kind(&self) -> &'static str {
        match self {
            Value::String(_) => "a string",
            Value::Integer(..) => "an integer",
            Value::Float(..) => "a float",
            Value::Boolean(_) => "a boolean",
            Value::Array(_) => "an array",
            Value::Table(_) => "a table",
        }
    }
}

/// A parse error, with its line (1-based).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub line: usize,
    pub message: String,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for Error {}

/// Parses a whole TOML document.
pub fn parse(text: &str) -> Result<Table, Error> {
    let mut p = Parser { s: text.as_bytes(), i: 0, line: 1 };
    let mut root = Table { defined: true, ..Table::default() };
    // the path of the current table, and whether it is an array of tables' last element
    let mut current: Vec<String> = Vec::new();
    loop {
        p.skip_ws_comments_newlines()?;
        if p.eof() {
            break;
        }
        if p.peek() == Some(b'[') {
            let array = p.peek_at(1) == Some(b'[');
            p.i += if array { 2 } else { 1 };
            p.skip_ws();
            let path = p.key_path()?;
            p.skip_ws();
            if !p.eat(b']') || (array && !p.eat(b']')) {
                return Err(p.err("expected ] after a table name"));
            }
            p.end_of_line()?;
            let line = p.line;
            if array {
                push_array_table(&mut root, &path).map_err(|m| Error { line, message: m })?;
            } else {
                define_table(&mut root, &path).map_err(|m| Error { line, message: m })?;
            }
            current = path;
        } else {
            let line = p.line;
            let key = p.key_path()?;
            p.skip_ws();
            if !p.eat(b'=') {
                return Err(p.err("expected = after a key"));
            }
            p.skip_ws();
            let v = p.value()?;
            p.end_of_line()?;
            let t = table_at(&mut root, &current).map_err(|m| Error { line, message: m })?;
            insert_dotted(t, &key, v).map_err(|m| Error { line, message: m })?;
        }
    }
    Ok(root)
}

/// The table at `path` (an array of tables resolves to its last element).
fn table_at<'a>(root: &'a mut Table, path: &[String]) -> Result<&'a mut Table, String> {
    let mut t = root;
    for k in path {
        let v = t.get_mut(k).ok_or_else(|| format!("no table {k}"))?;
        t = match v {
            Value::Table(t) => t,
            Value::Array(a) => match a.last_mut() {
                Some(Value::Table(t)) => t,
                _ => return Err(format!("{k} is not a table")),
            },
            _ => return Err(format!("{k} is not a table")),
        };
    }
    Ok(t)
}

/// The table at `path`, implied tables created on the way.
fn walk<'a>(root: &'a mut Table, path: &[String]) -> Result<&'a mut Table, String> {
    let mut t = root;
    for (n, k) in path.iter().enumerate() {
        if t.get(k).is_none() {
            t.entries.push((k.clone(), Value::Table(Table::default())));
        }
        let v = t.get_mut(k).unwrap();
        t = match v {
            Value::Table(t) => {
                if t.sealed {
                    return Err(format!("{} is an inline table; it cannot be extended", path[..=n].join(".")));
                }
                t
            }
            Value::Array(a) => match a.last_mut() {
                Some(Value::Table(t)) => t,
                _ => return Err(format!("{} is an array, not a table", path[..=n].join("."))),
            },
            other => return Err(format!("{} is {}, not a table", path[..=n].join("."), other.kind())),
        };
    }
    Ok(t)
}

fn define_table(root: &mut Table, path: &[String]) -> Result<(), String> {
    let (last, parent) = path.split_last().expect("a table name");
    let p = walk(root, parent)?;
    match p.get_mut(last) {
        None => {
            p.entries.push((last.clone(), Value::Table(Table { defined: true, ..Table::default() })));
            Ok(())
        }
        Some(Value::Table(t)) if !t.defined && !t.sealed => {
            t.defined = true;
            Ok(())
        }
        Some(_) => Err(format!("table {} defined twice", path.join("."))),
    }
}

fn push_array_table(root: &mut Table, path: &[String]) -> Result<(), String> {
    let (last, parent) = path.split_last().expect("a table name");
    let p = walk(root, parent)?;
    let fresh = Value::Table(Table { defined: true, ..Table::default() });
    match p.get_mut(last) {
        None => {
            p.entries.push((last.clone(), Value::Array(vec![fresh])));
            Ok(())
        }
        Some(Value::Array(a)) if a.iter().all(|v| matches!(v, Value::Table(t) if !t.sealed)) && !a.is_empty() => {
            a.push(fresh);
            Ok(())
        }
        Some(_) => Err(format!("{} is not an array of tables", path.join("."))),
    }
}

fn insert_dotted(t: &mut Table, key: &[String], v: Value) -> Result<(), String> {
    let (last, parent) = key.split_last().expect("a key");
    let t = walk(t, parent)?;
    if t.get(last).is_some() {
        return Err(format!("key {} defined twice", key.join(".")));
    }
    t.entries.push((last.clone(), v));
    Ok(())
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
    line: usize,
}

impl Parser<'_> {
    fn err(&self, m: &str) -> Error {
        Error { line: self.line, message: m.to_string() }
    }

    fn eof(&self) -> bool {
        self.i >= self.s.len()
    }

    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }

    fn peek_at(&self, k: usize) -> Option<u8> {
        self.s.get(self.i + k).copied()
    }

    fn eat(&mut self, b: u8) -> bool {
        if self.peek() == Some(b) {
            self.i += 1;
            true
        } else {
            false
        }
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t')) {
            self.i += 1;
        }
    }

    fn skip_comment(&mut self) {
        if self.peek() == Some(b'#') {
            while !matches!(self.peek(), None | Some(b'\n')) {
                self.i += 1;
            }
        }
    }

    fn newline(&mut self) -> bool {
        if self.peek() == Some(b'\n') {
            self.i += 1;
            self.line += 1;
            true
        } else if self.peek() == Some(b'\r') && self.peek_at(1) == Some(b'\n') {
            self.i += 2;
            self.line += 1;
            true
        } else {
            false
        }
    }

    fn skip_ws_comments_newlines(&mut self) -> Result<(), Error> {
        loop {
            self.skip_ws();
            self.skip_comment();
            if !self.newline() {
                return Ok(());
            }
        }
    }

    /// The rest of a line after a value or header: blanks, a comment, then a newline or the end.
    fn end_of_line(&mut self) -> Result<(), Error> {
        self.skip_ws();
        self.skip_comment();
        if self.eof() || self.newline() {
            Ok(())
        } else {
            Err(self.err("unexpected text after a value"))
        }
    }

    fn key_path(&mut self) -> Result<Vec<String>, Error> {
        let mut path = vec![self.simple_key()?];
        loop {
            self.skip_ws();
            if self.peek() == Some(b'.') {
                self.i += 1;
                self.skip_ws();
                path.push(self.simple_key()?);
            } else {
                return Ok(path);
            }
        }
    }

    fn simple_key(&mut self) -> Result<String, Error> {
        match self.peek() {
            Some(b'"') => self.basic_string(),
            Some(b'\'') => self.literal_string(),
            _ => {
                let start = self.i;
                while matches!(self.peek(), Some(b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' | b'-')) {
                    self.i += 1;
                }
                if self.i == start {
                    return Err(self.err("expected a key"));
                }
                Ok(String::from_utf8_lossy(&self.s[start..self.i]).into_owned())
            }
        }
    }

    fn value(&mut self) -> Result<Value, Error> {
        match self.peek() {
            Some(b'"') => {
                if self.s[self.i..].starts_with(b"\"\"\"") {
                    self.multiline_basic().map(Value::String)
                } else {
                    self.basic_string().map(Value::String)
                }
            }
            Some(b'\'') => {
                if self.s[self.i..].starts_with(b"'''") {
                    self.multiline_literal().map(Value::String)
                } else {
                    self.literal_string().map(Value::String)
                }
            }
            Some(b'[') => self.array(),
            Some(b'{') => self.inline_table(),
            Some(b't') if self.s[self.i..].starts_with(b"true") => {
                self.i += 4;
                Ok(Value::Boolean(true))
            }
            Some(b'f') if self.s[self.i..].starts_with(b"false") => {
                self.i += 5;
                Ok(Value::Boolean(false))
            }
            Some(_) => self.number(),
            None => Err(self.err("expected a value")),
        }
    }

    fn array(&mut self) -> Result<Value, Error> {
        self.i += 1; // [
        let mut items = Vec::new();
        loop {
            self.skip_ws_comments_newlines()?;
            if self.eat(b']') {
                return Ok(Value::Array(items));
            }
            items.push(self.value()?);
            self.skip_ws_comments_newlines()?;
            if self.eat(b',') {
                continue;
            }
            self.skip_ws_comments_newlines()?;
            if self.eat(b']') {
                return Ok(Value::Array(items));
            }
            return Err(self.err("expected , or ] in an array"));
        }
    }

    fn inline_table(&mut self) -> Result<Value, Error> {
        self.i += 1; // {
        let mut t = Table { defined: true, ..Table::default() };
        self.skip_ws();
        if self.eat(b'}') {
            t.sealed = true;
            return Ok(Value::Table(t));
        }
        loop {
            self.skip_ws();
            let key = self.key_path()?;
            self.skip_ws();
            if !self.eat(b'=') {
                return Err(self.err("expected = in an inline table"));
            }
            self.skip_ws();
            let v = self.value()?;
            let line = self.line;
            insert_dotted(&mut t, &key, v).map_err(|m| Error { line, message: m })?;
            self.skip_ws();
            if self.eat(b',') {
                continue;
            }
            if self.eat(b'}') {
                seal(&mut t);
                return Ok(Value::Table(t));
            }
            return Err(self.err("expected , or } in an inline table"));
        }
    }

    fn escape(&mut self, out: &mut Vec<u8>) -> Result<(), Error> {
        let c = self.peek().ok_or_else(|| self.err("unterminated escape"))?;
        self.i += 1;
        match c {
            b'b' => out.push(8),
            b't' => out.push(b'\t'),
            b'n' => out.push(b'\n'),
            b'f' => out.push(12),
            b'r' => out.push(b'\r'),
            b'"' => out.push(b'"'),
            b'\\' => out.push(b'\\'),
            b'u' | b'U' => {
                let n = if c == b'u' { 4 } else { 8 };
                let hex = self.s.get(self.i..self.i + n).ok_or_else(|| self.err("short unicode escape"))?;
                let code = u32::from_str_radix(std::str::from_utf8(hex).unwrap_or("x"), 16).map_err(|_| self.err("bad unicode escape"))?;
                let ch = char::from_u32(code).ok_or_else(|| self.err("bad unicode scalar"))?;
                let mut buf = [0u8; 4];
                out.extend_from_slice(ch.encode_utf8(&mut buf).as_bytes());
                self.i += n;
            }
            _ => return Err(self.err("unknown escape")),
        }
        Ok(())
    }

    fn utf8(&self, v: Vec<u8>) -> Result<String, Error> {
        String::from_utf8(v).map_err(|_| self.err("invalid UTF-8"))
    }

    fn basic_string(&mut self) -> Result<String, Error> {
        self.i += 1; // "
        let mut out = Vec::new();
        loop {
            match self.peek() {
                None | Some(b'\n') => return Err(self.err("unterminated string")),
                Some(b'"') => {
                    self.i += 1;
                    return self.utf8(out);
                }
                Some(b'\\') => {
                    self.i += 1;
                    self.escape(&mut out)?;
                }
                Some(c) => {
                    out.push(c);
                    self.i += 1;
                }
            }
        }
    }

    fn literal_string(&mut self) -> Result<String, Error> {
        self.i += 1; // '
        let start = self.i;
        while !matches!(self.peek(), None | Some(b'\'' | b'\n')) {
            self.i += 1;
        }
        if !self.eat(b'\'') {
            return Err(self.err("unterminated literal string"));
        }
        self.utf8(self.s[start..self.i - 1].to_vec())
    }

    fn multiline_basic(&mut self) -> Result<String, Error> {
        self.i += 3;
        self.newline(); // a newline right after the delimiter is trimmed
        let mut out = Vec::new();
        loop {
            if self.s[self.i..].starts_with(b"\"\"\"") {
                // up to two quotes may end the content
                let mut q = 3;
                while self.peek_at(q) == Some(b'"') && q < 5 {
                    q += 1;
                }
                out.extend(std::iter::repeat_n(b'"', q - 3));
                self.i += q;
                return self.utf8(out);
            }
            match self.peek() {
                None => return Err(self.err("unterminated multi-line string")),
                Some(b'\\') => {
                    self.i += 1;
                    // a line-ending backslash trims the newline and the next blanks
                    let save = self.i;
                    self.skip_ws();
                    if self.newline() {
                        loop {
                            self.skip_ws();
                            if !self.newline() {
                                break;
                            }
                        }
                    } else {
                        self.i = save;
                        self.escape(&mut out)?;
                    }
                }
                Some(b'\n') => {
                    self.newline();
                    out.push(b'\n');
                }
                Some(c) => {
                    out.push(c);
                    self.i += 1;
                }
            }
        }
    }

    fn multiline_literal(&mut self) -> Result<String, Error> {
        self.i += 3;
        self.newline();
        let mut out = Vec::new();
        loop {
            if self.s[self.i..].starts_with(b"'''") {
                let mut q = 3;
                while self.peek_at(q) == Some(b'\'') && q < 5 {
                    q += 1;
                }
                out.extend(std::iter::repeat_n(b'\'', q - 3));
                self.i += q;
                return self.utf8(out);
            }
            match self.peek() {
                None => return Err(self.err("unterminated multi-line string")),
                Some(b'\n') => {
                    self.newline();
                    out.push(b'\n');
                }
                Some(c) => {
                    out.push(c);
                    self.i += 1;
                }
            }
        }
    }

    fn number(&mut self) -> Result<Value, Error> {
        let start = self.i;
        while matches!(self.peek(), Some(b'0'..=b'9' | b'a'..=b'z' | b'A'..=b'Z' | b'+' | b'-' | b'.' | b'_' | b':')) {
            self.i += 1;
        }
        let raw = std::str::from_utf8(&self.s[start..self.i]).unwrap_or("");
        if raw.is_empty() {
            return Err(self.err("expected a value"));
        }
        if raw.contains(':') || (raw.len() >= 10 && raw.as_bytes()[4] == b'-' && raw.as_bytes()[..4].iter().all(u8::is_ascii_digit)) {
            return Err(self.err("dates and times are not supported"));
        }
        let bad = || self.err(&format!("not a number: {raw}"));
        let (sign, body) = match raw.as_bytes()[0] {
            b'+' => (1i64, &raw[1..]),
            b'-' => (-1i64, &raw[1..]),
            _ => (1i64, raw),
        };
        match body {
            "inf" => return Ok(Value::Float(sign as f64 * f64::INFINITY, raw.into())),
            "nan" => return Ok(Value::Float(f64::NAN, raw.into())),
            _ => {}
        }
        // underscores only between digits
        let b = body.as_bytes();
        for (k, &c) in b.iter().enumerate() {
            if c == b'_' && !(k > 0 && k + 1 < b.len() && b[k - 1].is_ascii_alphanumeric() && b[k + 1].is_ascii_alphanumeric()) {
                return Err(bad());
            }
        }
        let clean: String = raw.chars().filter(|&c| c != '_').collect();
        let digits: String = body.chars().filter(|&c| c != '_').collect();
        for (prefix, radix) in [("0x", 16), ("0o", 8), ("0b", 2)] {
            if let Some(d) = digits.strip_prefix(prefix) {
                if sign != 1 || raw.starts_with('+') {
                    return Err(bad());
                }
                let v = i64::from_str_radix(d, radix).map_err(|_| bad())?;
                return Ok(Value::Integer(v, clean));
            }
        }
        let is_float = digits.contains(['.', 'e', 'E']);
        // no leading zeros in the integer part
        let int_part: &str = digits.split(['.', 'e', 'E']).next().unwrap_or("");
        if int_part.is_empty() || !int_part.bytes().all(|c| c.is_ascii_digit()) || (int_part.len() > 1 && int_part.starts_with('0')) {
            return Err(bad());
        }
        if is_float {
            // a dot must have digits on both sides
            if let Some(k) = digits.find('.') {
                let after = digits.as_bytes().get(k + 1);
                if !after.is_some_and(u8::is_ascii_digit) {
                    return Err(bad());
                }
            }
            let v: f64 = clean.parse().map_err(|_| bad())?;
            Ok(Value::Float(v, clean))
        } else {
            let v: i64 = clean.parse().map_err(|_| bad())?;
            Ok(Value::Integer(v, clean))
        }
    }
}

/// Seals an inline table and the tables inside it.
fn seal(t: &mut Table) {
    t.sealed = true;
    for (_, v) in &mut t.entries {
        if let Value::Table(s) = v {
            seal(s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_and_tables() {
        let t = parse(
            "# comment\npreset = \"workflow\" # trailing\n\n[cluster]\nresolution = [0.8, 1.2,\n  2.0, # last\n]\nseed = 12_345\nkeep = true\nhex = 0xff\n\n[neighbors]\nprune = '1/15'\nmetric.name = \"euclidean\"\nx = -1.5e-3\n",
        )
        .unwrap();
        assert_eq!(t.get("preset").unwrap().as_str(), Some("workflow"));
        let c = t.get("cluster").unwrap().as_table().unwrap();
        let r = c.get("resolution").unwrap().as_array().unwrap();
        assert_eq!(r.len(), 3);
        assert_eq!(r[1].number_text(), Some("1.2"));
        assert_eq!(r[2].as_f64(), Some(2.0));
        assert_eq!(c.get("seed").unwrap().as_i64(), Some(12345));
        assert_eq!(c.get("seed").unwrap().number_text(), Some("12345"));
        assert_eq!(c.get("keep").unwrap().as_bool(), Some(true));
        assert_eq!(c.get("hex").unwrap().as_i64(), Some(255));
        assert_eq!(t.path("neighbors.prune").unwrap().as_str(), Some("1/15"));
        assert_eq!(t.path("neighbors.metric.name").unwrap().as_str(), Some("euclidean"));
        assert_eq!(t.path("neighbors.x").unwrap().as_f64(), Some(-1.5e-3));
        assert_eq!(c.keys().collect::<Vec<_>>(), vec!["resolution", "seed", "keep", "hex"]);
    }

    #[test]
    fn strings_and_inline_tables() {
        let t = parse("a = \"tab\\tq\\\"u\\u00e9\"\nb = '''\nraw \\n\n'''\nc = \"\"\"\nx \\\n   y\"\"\"\nd = { k = 1, s.t = \"v\" }\ne = [ { x = 1 }, { x = 2 } ]\n").unwrap();
        assert_eq!(t.get("a").unwrap().as_str(), Some("tab\tq\"u\u{e9}"));
        assert_eq!(t.get("b").unwrap().as_str(), Some("raw \\n\n"));
        assert_eq!(t.get("c").unwrap().as_str(), Some("x y"));
        assert_eq!(t.path("d.s.t").unwrap().as_str(), Some("v"));
        assert_eq!(t.get("e").unwrap().as_array().unwrap()[1].as_table().unwrap().get("x").unwrap().as_i64(), Some(2));
    }

    #[test]
    fn arrays_of_tables() {
        let t = parse("[[file]]\npath = \"a\"\n[[file]]\npath = \"b\"\n[file.sub]\nx = 1\n[[run]]\nid = 1\n").unwrap();
        let f = t.get("file").unwrap().as_array().unwrap();
        assert_eq!(f.len(), 2);
        assert_eq!(f[1].as_table().unwrap().path("sub.x").unwrap().as_i64(), Some(1));
        assert_eq!(t.get("run").unwrap().as_array().unwrap().len(), 1);
    }

    #[test]
    fn errors() {
        for bad in [
            "a = 1\na = 2\n",
            "[t]\nx = 1\n[t]\ny = 2\n",
            "a = { x = 1 }\n[a]\ny = 2\n",
            "a = 01\n",
            "a = 1.\n",
            "a = 1979-05-27\n",
            "a = \"open\n",
            "a = [1, 2\n",
            "a = 1 b\n",
            "= 1\n",
            "a.b = 1\na = 2\n",
            "x = 1__0\n",
        ] {
            assert!(parse(bad).is_err(), "accepted: {bad:?}");
        }
        let e = parse("a = 1\n\nb = 1 x\n").unwrap_err();
        assert_eq!(e.line, 3);
        // a table implied by a sub-table may be defined once later
        assert!(parse("[a.b]\nx = 1\n[a]\ny = 2\n").is_ok());
    }

    #[test]
    fn special_floats() {
        let t = parse("a = inf\nb = -inf\nc = nan\nd = +1e2\n").unwrap();
        assert_eq!(t.get("a").unwrap().as_f64(), Some(f64::INFINITY));
        assert_eq!(t.get("b").unwrap().as_f64(), Some(f64::NEG_INFINITY));
        assert!(t.get("c").unwrap().as_f64().unwrap().is_nan());
        assert_eq!(t.get("d").unwrap().as_f64(), Some(100.0));
    }

    #[test]
    fn the_repository_s_tolerances() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../validation/tolerances.toml");
        let t = parse(&std::fs::read_to_string(path).unwrap()).unwrap();
        assert!(t.path("doubletfinder.floor.kernel.S3_sdev.max_rel").is_some());
        assert!(t.path("clustering.floor.kernel").is_some());
    }
}

//! Small helpers shared by the subcommands: `--name value` arguments, run records (time, peak
//! memory) and errors.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// `--name value` arguments.
pub struct Args {
    values: BTreeMap<String, String>,
}

impl Args {
    pub fn parse(args: &[String]) -> Result<Args, String> {
        let mut values = BTreeMap::new();
        let mut i = 0;
        while i < args.len() {
            let a = &args[i];
            let Some(name) = a.strip_prefix("--") else {
                return Err(format!("unexpected argument '{a}'"));
            };
            let v = args.get(i + 1).ok_or_else(|| format!("--{name} needs a value"))?;
            if values.insert(name.to_string(), v.clone()).is_some() {
                return Err(format!("--{name} given twice"));
            }
            i += 2;
        }
        Ok(Args { values })
    }

    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str)
    }

    pub fn required(&self, name: &str) -> Result<&str, String> {
        self.get(name).ok_or_else(|| format!("--{name} is required"))
    }

    pub fn path(&self, name: &str) -> Result<PathBuf, String> {
        self.required(name).map(PathBuf::from)
    }

    pub fn number<T: std::str::FromStr>(&self, name: &str, default: T) -> Result<T, String> {
        match self.get(name) {
            None => Ok(default),
            Some(v) => v.parse().map_err(|_| format!("--{name}: invalid value '{v}'")),
        }
    }

    /// Fails on any argument not in `known`.
    pub fn only(&self, known: &[&str]) -> Result<(), String> {
        for k in self.values.keys() {
            if !known.contains(&k.as_str()) {
                return Err(format!("unknown option --{k}"));
            }
        }
        Ok(())
    }
}

/// The number of threads: `--threads`, else `SLURM_CPUS_PER_TASK`, else 1.
pub fn threads(args: &Args) -> Result<usize, String> {
    let default = std::env::var("SLURM_CPUS_PER_TASK").ok().and_then(|v| v.parse().ok()).unwrap_or(1);
    let t: usize = args.number("threads", default)?;
    Ok(t.max(1))
}

/// Peak resident set size of this process so far (`VmHWM`), in KiB.
pub fn peak_rss_kib() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines().find_map(|l| l.strip_prefix("VmHWM:")).and_then(|v| v.trim().trim_end_matches("kB").trim().parse().ok())
}

/// The current time in UTC, ISO 8601.
pub fn utc_now() -> String {
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let (days, rem) = ((secs / 86400) as i64, secs % 86400);
    // civil date from days since 1970-01-01 (Howard Hinnant's algorithm)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + (m <= 2) as i64;
    format!("{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z", rem / 3600, (rem / 60) % 60, rem % 60)
}

/// Elapsed time and peak memory after each phase of a run.
pub struct Phases {
    start: Instant,
    pub marks: Vec<(String, f64, Option<u64>)>,
}

impl Phases {
    pub fn new() -> Phases {
        Phases { start: Instant::now(), marks: Vec::new() }
    }

    pub fn mark(&mut self, name: &str) {
        self.marks.push((name.to_string(), self.start.elapsed().as_secs_f64(), peak_rss_kib()));
    }

    pub fn elapsed(&self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }

    pub fn to_json(&self) -> scio::Json {
        scio::Json::Object(
            self.marks
                .iter()
                .map(|(n, t, r)| {
                    (
                        n.clone(),
                        scio::Json::obj(vec![
                            ("elapsed_s", scio::Json::Num(*t)),
                            ("peak_rss_kib", r.map(|v| scio::Json::Int(v as i64)).unwrap_or(scio::Json::Null)),
                        ]),
                    )
                })
                .collect(),
        )
    }
}

/// An environment variable as JSON (null when unset).
pub fn env_json(name: &str) -> scio::Json {
    std::env::var(name).map(scio::Json::Str).unwrap_or(scio::Json::Null)
}

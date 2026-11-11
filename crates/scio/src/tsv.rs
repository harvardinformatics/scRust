//! Tab-separated tables with a header line, as the reference scripts write them (`tsv_write`
//! in `validation/R/lib/io.R`): no quoting, `\n` line ends, one header line.

use std::fs;
use std::io::{self, Write};
use std::path::Path;

/// A table: column names and rows of fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl Table {
    pub fn read(path: &Path) -> io::Result<Table> {
        let text = fs::read_to_string(path)
            .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
        Table::parse(&text).ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, format!("{}: empty table", path.display())))
    }

    pub fn parse(text: &str) -> Option<Table> {
        let mut lines = text.lines();
        let header = lines.next()?.split('\t').map(str::to_owned).collect();
        let rows = lines.map(|l| l.split('\t').map(str::to_owned).collect()).collect();
        Some(Table { header, rows })
    }

    /// The index of a column.
    pub fn col(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|h| h == name)
    }

    /// The table as text: header and rows joined by tabs, each line ending in `\n`.
    pub fn to_text(&self) -> String {
        let mut s = String::new();
        for line in std::iter::once(&self.header).chain(&self.rows) {
            s.push_str(&line.join("\t"));
            s.push('\n');
        }
        s
    }

    pub fn write(&self, path: &Path) -> io::Result<()> {
        fs::File::create(path)?.write_all(self.to_text().as_bytes())
    }
}

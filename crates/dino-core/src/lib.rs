//! Shared suite metadata for Dino Tools.

/// How a tool currently relates to the suite workspace.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ToolStatus {
    /// The tool exists outside this workspace and should stay independently developed.
    External,
    /// The tool name and role are reserved, but no suite-owned implementation exists yet.
    Planned,
    /// The tool is implemented inside this workspace.
    Workspace,
}

impl ToolStatus {
    /// Stable display label used by the CLI and docs.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::External => "external",
            Self::Planned => "planned",
            Self::Workspace => "workspace",
        }
    }
}

/// Where the implementation for a Dino Tools tool currently lives.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ToolSource {
    /// The tool is still developed in a standalone repository.
    ExternalRepo {
        local_path: &'static str,
        repository: &'static str,
    },
    /// The tool is implemented as part of this workspace.
    WorkspaceCrate { crate_name: &'static str },
    /// The tool name is reserved, but there is no implementation yet.
    Planned,
}

impl ToolSource {
    #[must_use]
    pub const fn kind(self) -> &'static str {
        match self {
            Self::ExternalRepo { .. } => "external repo",
            Self::WorkspaceCrate { .. } => "workspace crate",
            Self::Planned => "planned",
        }
    }
}

/// A user-facing Dino Tools tool entry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ToolDescriptor {
    pub name: &'static str,
    pub role: &'static str,
    pub status: ToolStatus,
    pub source: ToolSource,
    pub command: Option<&'static str>,
    pub suite_contracts: &'static [&'static str],
}

impl ToolDescriptor {
    #[must_use]
    pub const fn external(
        name: &'static str,
        role: &'static str,
        local_path: &'static str,
        repository: &'static str,
        command: &'static str,
        suite_contracts: &'static [&'static str],
    ) -> Self {
        Self {
            name,
            role,
            status: ToolStatus::External,
            source: ToolSource::ExternalRepo {
                local_path,
                repository,
            },
            command: Some(command),
            suite_contracts,
        }
    }

    #[must_use]
    pub const fn workspace(
        name: &'static str,
        role: &'static str,
        crate_name: &'static str,
        command: &'static str,
        suite_contracts: &'static [&'static str],
    ) -> Self {
        Self {
            name,
            role,
            status: ToolStatus::Workspace,
            source: ToolSource::WorkspaceCrate { crate_name },
            command: Some(command),
            suite_contracts,
        }
    }

    #[must_use]
    pub const fn planned(name: &'static str, role: &'static str) -> Self {
        Self {
            name,
            role,
            status: ToolStatus::Planned,
            source: ToolSource::Planned,
            command: None,
            suite_contracts: &[],
        }
    }
}

/// The initial tool registry. Existing tools are tracked as external on purpose.
pub const TOOL_REGISTRY: &[ToolDescriptor] = &[
    ToolDescriptor::external(
        "trex",
        "assembler / heavy de novo engine",
        "/home/jake/Projects/Trex",
        "local",
        "trex",
        &["read ingest boundary", "assembly pipeline caller"],
    ),
    ToolDescriptor::workspace(
        "dino-seq",
        "FASTQ/FASTA parsing and ingest",
        "dino-seq",
        "dino-seq",
        &[
            "raw/gzip/BGZF FASTQ opener",
            "raw/gzip/BGZF FASTA opener",
            "streaming batch readers",
            "FASTA index and reference chunk helpers",
        ],
    ),
    ToolDescriptor::planned("velociraptor", "fast search, sketching, and prefiltering"),
    ToolDescriptor::planned("ankylosaurus", "QC, validation, and defensive filtering"),
    ToolDescriptor::planned("stegosaurus", "graph layout and scaffolding"),
    ToolDescriptor::planned("triceratops", "comparison, reconciliation, and consensus"),
    ToolDescriptor::planned("brachiosaurus", "long-read and large-reference workflows"),
    ToolDescriptor::planned("archaeopteryx", "import/export and compatibility bridges"),
];

/// Find a registered tool by its command name.
#[must_use]
pub fn find_tool(name: &str) -> Option<&'static ToolDescriptor> {
    TOOL_REGISTRY.iter().find(|tool| tool.name == name)
}

#[cfg(test)]
mod tests {
    use super::{find_tool, ToolSource, ToolStatus, TOOL_REGISTRY};

    #[test]
    fn registry_tracks_trex_as_external_and_dino_seq_as_workspace() {
        assert_eq!(
            find_tool("trex").map(|tool| tool.status),
            Some(ToolStatus::External)
        );
        assert_eq!(
            find_tool("dino-seq").map(|tool| tool.status),
            Some(ToolStatus::Workspace)
        );
    }

    #[test]
    fn registry_names_are_unique() {
        for (index, left) in TOOL_REGISTRY.iter().enumerate() {
            for right in &TOOL_REGISTRY[index + 1..] {
                assert_ne!(left.name, right.name);
            }
        }
    }

    #[test]
    fn dino_seq_records_workspace_contracts() {
        let tool = find_tool("dino-seq");

        assert_eq!(
            tool.map(|descriptor| descriptor.command),
            Some(Some("dino-seq"))
        );
        assert_eq!(
            tool.map(|descriptor| descriptor.source),
            Some(ToolSource::WorkspaceCrate {
                crate_name: "dino-seq",
            })
        );
        assert_eq!(
            tool.map(|descriptor| descriptor
                .suite_contracts
                .contains(&"streaming batch readers")),
            Some(true)
        );
    }
}

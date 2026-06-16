use std::env;
use std::process::ExitCode;

use dino_core::{find_tool, ToolSource, TOOL_REGISTRY};

fn main() -> ExitCode {
    let mut args = env::args().skip(1);
    match args.next().as_deref() {
        None | Some("-h" | "--help" | "help") => {
            print_help();
            ExitCode::SUCCESS
        }
        Some("list" | "status") => {
            print_tools();
            ExitCode::SUCCESS
        }
        Some(tool_name) => match find_tool(tool_name) {
            Some(tool) => {
                println!(
                    "{} is currently {}: {}",
                    tool.name,
                    tool.status.as_str(),
                    tool.role
                );
                print_tool_source(tool.source);
                if let Some(command) = tool.command {
                    println!("command: {command}");
                }
                if !tool.suite_contracts.is_empty() {
                    println!("suite contracts:");
                    for contract in tool.suite_contracts {
                        println!("  - {contract}");
                    }
                }
                println!("Standalone tool dispatch is intentionally not wired yet.");
                ExitCode::SUCCESS
            }
            None => {
                eprintln!("unknown dino command or tool: {tool_name}");
                eprintln!("run `dino list` to see registered names");
                ExitCode::FAILURE
            }
        },
    }
}

fn print_help() {
    println!("Dino Tools bioinformatics suite");
    println!();
    println!("Usage:");
    println!("  dino list");
    println!("  dino status");
    println!("  dino <tool>");
    println!();
    println!("Existing tools stay standalone until explicitly promoted.");
}

fn print_tools() {
    println!("{:<15} {:<10} {:<13} role", "tool", "status", "source");
    println!("{:<15} {:<10} {:<13} ----", "----", "------", "------");
    for tool in TOOL_REGISTRY {
        println!(
            "{:<15} {:<10} {:<13} {}",
            tool.name,
            tool.status.as_str(),
            tool.source.kind(),
            tool.role
        );
    }
}

fn print_tool_source(source: ToolSource) {
    match source {
        ToolSource::ExternalRepo {
            local_path,
            repository,
        } => {
            println!("source: external repo");
            println!("local path: {local_path}");
            println!("repository: {repository}");
        }
        ToolSource::WorkspaceCrate { crate_name } => {
            println!("source: workspace crate");
            println!("crate: {crate_name}");
        }
        ToolSource::Planned => {
            println!("source: planned");
        }
    }
}

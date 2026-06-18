use std::error::Error;
use std::ffi::OsStr;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_dir(name: &str) -> Result<PathBuf, Box<dyn Error>> {
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let dir =
        std::env::temp_dir().join(format!("dino-quant-{name}-{}-{nanos}", std::process::id()));
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

fn run_dino_quant<I, S>(args: I) -> Result<Output, Box<dyn Error>>
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    Ok(Command::new(env!("CARGO_BIN_EXE_dino-quant"))
        .args(args)
        .output()?)
}

fn assert_success(output: &Output) -> Result<(), Box<dyn Error>> {
    if output.status.success() {
        return Ok(());
    }
    Err(format!(
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
    .into())
}

#[test]
fn minimizer_cache_rejects_same_path_reference_content_change() -> Result<(), Box<dyn Error>> {
    let dir = temp_dir("cache-fingerprint")?;
    let ref_path = dir.join("ref.fa");
    let reads_path = dir.join("reads.fastq");
    let cache_path = dir.join("ref.dqmi");
    fs::write(&ref_path, b">ref\nACGTACGTGGGGTTTT\n")?;
    fs::write(&reads_path, b"@r1\nACGTACGT\n+\nFFFFFFFF\n")?;

    let args = [
        OsStr::new("emit-candidates"),
        ref_path.as_os_str(),
        reads_path.as_os_str(),
        OsStr::new("--retrieval"),
        OsStr::new("minimizer"),
        OsStr::new("--top-k"),
        OsStr::new("2"),
        OsStr::new("--window"),
        OsStr::new("8"),
        OsStr::new("--stride"),
        OsStr::new("4"),
        OsStr::new("--k"),
        OsStr::new("3"),
        OsStr::new("--dim"),
        OsStr::new("16"),
        OsStr::new("--candidate-limit"),
        OsStr::new("8"),
        OsStr::new("--minimizer-k"),
        OsStr::new("3"),
        OsStr::new("--minimizer-window"),
        OsStr::new("2"),
        OsStr::new("--reference-cache"),
        cache_path.as_os_str(),
    ];

    let first = run_dino_quant(args)?;
    assert_success(&first)?;
    assert!(String::from_utf8_lossy(&first.stderr).contains("cache_hit=false"));

    let second = run_dino_quant(args)?;
    assert_success(&second)?;
    assert!(String::from_utf8_lossy(&second.stderr).contains("cache_hit=true"));

    fs::write(&ref_path, b">ref\nTTTTCCCCAAAAGGGG\n")?;
    let changed = run_dino_quant(args)?;
    assert_success(&changed)?;
    assert!(String::from_utf8_lossy(&changed.stderr).contains("cache_hit=false"));

    fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn candidate_reference_uses_target_id_for_duplicate_names() -> Result<(), Box<dyn Error>> {
    let dir = temp_dir("target-id")?;
    let ref_path = dir.join("ref.fa");
    let candidates_path = dir.join("candidates.tsv");
    fs::write(
        &ref_path,
        b">dup first copy\nAAAAAAAA\n>dup second copy\nCCCCCCCC\n",
    )?;
    fs::write(
        &candidates_path,
        b"query_name\tquery_len\trank\ttarget_name\ttarget_start\ttarget_end\tlinear_start\tlinear_end\tscore\ttarget_id\nq\t8\t1\tdup\t2\t6\t10\t14\t1.000000\t1\n",
    )?;

    let output = run_dino_quant([
        OsStr::new("emit-candidate-reference"),
        ref_path.as_os_str(),
        candidates_path.as_os_str(),
        OsStr::new("--mask-reference"),
        OsStr::new("--padding"),
        OsStr::new("0"),
    ])?;
    assert_success(&output)?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout, ">dup\nNNNNNNNN\n>dup\nNNCCCCNN\n");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("records=2"));
    assert!(stderr.contains("bases=16"));

    fs::remove_dir_all(&dir)?;
    Ok(())
}

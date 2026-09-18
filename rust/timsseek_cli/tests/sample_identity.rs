use std::process::Command;

#[test]
fn cli_rejects_duplicate_samples_before_opening_inputs_with_or_without_overwrite() {
    for overwrite in [false, true] {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().join("output");
        let mut command = Command::new(env!("CARGO_BIN_EXE_timsseek"));
        command.args([
            "--raw-inputs",
            "s3://bucket/run.d",
            "--raw-inputs",
            "s3://bucket/run.d.idx/",
            "--speclib-uri",
            "absent-library.mzspeclib.txt",
            "--output-uri",
            output.to_str().unwrap(),
            "--log-path",
            "-",
        ]);
        if overwrite {
            command.arg("--overwrite");
        }
        let result = command.output().unwrap();
        assert_eq!(result.status.code(), Some(1));
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(stderr.contains("Duplicate sample_id"), "{stderr}");
        assert!(stderr.contains("s3://bucket/run.d"), "{stderr}");
        assert!(stderr.contains("s3://bucket/run.d.idx/"), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        assert!(!output.exists());
    }
}

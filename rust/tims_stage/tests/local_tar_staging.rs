use std::io::Write;
use tempfile::TempDir;
use tims_stage::{
    PerRunTempdir,
    SourceSpec,
    StagingBackend,
    StagingConfig,
};

fn write_ustar_header(buf: &mut Vec<u8>, name: &str, size: u64, typeflag: u8) {
    let start = buf.len();
    buf.resize(start + 512, 0);
    let blk = &mut buf[start..start + 512];
    let nb = name.as_bytes();
    blk[..nb.len()].copy_from_slice(nb);
    let s = format!("{:011o}", size);
    blk[124..124 + s.len()].copy_from_slice(s.as_bytes());
    blk[156] = typeflag;
    for b in &mut blk[148..156] {
        *b = b' ';
    }
}

fn pad(buf: &mut Vec<u8>) {
    let p = (512 - buf.len() % 512) % 512;
    buf.extend(std::iter::repeat(0u8).take(p));
}

#[test]
fn stage_local_tar_materializes_sample_d_wrapper() {
    let dir = TempDir::new().unwrap();
    let tar_path = dir.path().join("sample.d.tar");
    let mut tar = Vec::new();
    let tdf_body = vec![1u8; 7];
    let tdf_bin_body = vec![2u8; 1024];
    write_ustar_header(&mut tar, "analysis.tdf", tdf_body.len() as u64, b'0');
    tar.extend_from_slice(&tdf_body);
    pad(&mut tar);
    write_ustar_header(
        &mut tar,
        "analysis.tdf_bin",
        tdf_bin_body.len() as u64,
        b'0',
    );
    tar.extend_from_slice(&tdf_bin_body);
    pad(&mut tar);
    tar.extend(std::iter::repeat(0u8).take(1024));
    std::fs::File::create(&tar_path)
        .unwrap()
        .write_all(&tar)
        .unwrap();

    let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
    let spec = SourceSpec::LocalTar { path: tar_path };
    let staged = backend.stage(&spec).unwrap();

    assert!(staged.as_ref().join("analysis.tdf").is_file());
    assert!(staged.as_ref().join("analysis.tdf_bin").is_file());
    assert_eq!(
        std::fs::read(staged.as_ref().join("analysis.tdf")).unwrap(),
        tdf_body
    );
    assert_eq!(
        std::fs::read(staged.as_ref().join("analysis.tdf_bin")).unwrap(),
        tdf_bin_body
    );
}

#[test]
fn stage_local_tar_tempdir_vanishes_on_drop() {
    let dir = TempDir::new().unwrap();
    let tar_path = dir.path().join("sample.d.tar");
    let mut tar = Vec::new();
    write_ustar_header(&mut tar, "analysis.tdf", 3, b'0');
    tar.extend_from_slice(&[1, 2, 3]);
    pad(&mut tar);
    write_ustar_header(&mut tar, "analysis.tdf_bin", 3, b'0');
    tar.extend_from_slice(&[4, 5, 6]);
    pad(&mut tar);
    tar.extend(std::iter::repeat(0u8).take(1024));
    std::fs::File::create(&tar_path)
        .unwrap()
        .write_all(&tar)
        .unwrap();

    let backend = PerRunTempdir::new(StagingConfig::default()).unwrap();
    let spec = SourceSpec::LocalTar { path: tar_path };
    let tempdir_path_copy;
    {
        let staged = backend.stage(&spec).unwrap();
        tempdir_path_copy = staged.as_ref().parent().unwrap().to_path_buf();
        assert!(tempdir_path_copy.is_dir());
    }
    assert!(
        !tempdir_path_copy.is_dir(),
        "tempdir should be gone after drop"
    );
}

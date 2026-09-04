use filetime::FileTime;
use std::time::{
    Duration,
    SystemTime,
};
use tempfile::TempDir;
use tims_stage::{
    PerRunTempdir,
    StagingConfig,
};

fn backdate(path: &std::path::Path, hours_ago: u64) {
    let ts = SystemTime::now() - Duration::from_secs(hours_ago * 3600);
    filetime::set_file_mtime(path, FileTime::from_system_time(ts)).unwrap();
}

#[test]
fn sweep_removes_old_unlocked_tempdirs() {
    let root = TempDir::new().unwrap();
    let stale = root.path().join("timsseek-staging-ABC123");
    std::fs::create_dir(&stale).unwrap();
    backdate(&stale, 48);
    let cfg = StagingConfig {
        tempdir_root: Some(root.path().to_path_buf()),
        stale_sweep_age_hours: 24,
        ..StagingConfig::default()
    };
    let _backend = PerRunTempdir::new(cfg).unwrap();
    assert!(
        !stale.is_dir(),
        "48h-old unlocked dir should have been swept"
    );
}

#[test]
fn sweep_preserves_recent_tempdirs() {
    let root = TempDir::new().unwrap();
    let fresh = root.path().join("timsseek-staging-FRESH");
    std::fs::create_dir(&fresh).unwrap();
    // mtime defaults to now; no backdating.
    let cfg = StagingConfig {
        tempdir_root: Some(root.path().to_path_buf()),
        stale_sweep_age_hours: 24,
        ..StagingConfig::default()
    };
    let _backend = PerRunTempdir::new(cfg).unwrap();
    assert!(fresh.is_dir());
}

#[test]
fn sweep_preserves_locked_tempdirs() {
    let root = TempDir::new().unwrap();
    let locked = root.path().join("timsseek-staging-LIVE");
    std::fs::create_dir(&locked).unwrap();
    std::fs::File::create(locked.join(".lock")).unwrap();
    backdate(&locked, 36);
    let cfg = StagingConfig {
        tempdir_root: Some(root.path().to_path_buf()),
        stale_sweep_age_hours: 24,
        ..StagingConfig::default()
    };
    let _backend = PerRunTempdir::new(cfg).unwrap();
    assert!(
        locked.is_dir(),
        "locked dir should survive sweep when mtime < 2x threshold"
    );
}

#[test]
fn sweep_ignores_unrelated_dirs() {
    let root = TempDir::new().unwrap();
    let unrelated = root.path().join("my-important-data");
    std::fs::create_dir(&unrelated).unwrap();
    backdate(&unrelated, 48);
    let cfg = StagingConfig {
        tempdir_root: Some(root.path().to_path_buf()),
        stale_sweep_age_hours: 24,
        ..StagingConfig::default()
    };
    let _backend = PerRunTempdir::new(cfg).unwrap();
    assert!(
        unrelated.is_dir(),
        "non-namespaced dir should be ignored by sweep"
    );
}

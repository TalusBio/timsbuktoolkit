use crate::common::transport_err;
use crate::error::StageError;
use crate::uri::{
    is_remote_uri,
    split_uri,
};
use std::path::Path;
use timscentroid::StorageProvider;

pub fn upload_file(local: &Path, dest_uri: &str) -> Result<(), StageError> {
    let bytes = std::fs::read(local).map_err(StageError::Io)?;
    if is_remote_uri(dest_uri) {
        let (loc, key) = split_uri(dest_uri)?;
        // Writes go through `new` (creates parent dir for local).
        let provider = StorageProvider::new(loc).map_err(transport_err(dest_uri))?;
        provider
            .write_bytes(&key, bytes)
            .map_err(transport_err(dest_uri))?;
    } else {
        if let Some(parent) = Path::new(dest_uri).parent() {
            std::fs::create_dir_all(parent).map_err(StageError::Io)?;
        }
        std::fs::write(dest_uri, &bytes).map_err(StageError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn copies_local_to_local() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("src.bin");
        let dst = dir.path().join("subdir/dst.bin");
        std::fs::write(&src, b"payload").unwrap();
        upload_file(&src, dst.to_str().unwrap()).unwrap();
        assert_eq!(std::fs::read(&dst).unwrap(), b"payload");
    }
}

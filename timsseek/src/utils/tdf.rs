use rusqlite::{
    Connection,
    Result,
};
use std::sync::Arc;

pub fn get_ms1_frame_times_ms(tdf_path: &str) -> Result<Arc<[u32]>> {
    let conn = Connection::open(tdf_path)?;
    let mut stmt = conn.prepare("SELECT Time FROM Frames WHERE MsMsType == 0")?;

    let times: Vec<f32> = stmt
        .query_map([], |row| row.get(0))?
        .collect::<Result<Vec<f32>>>()?;
    let times: Vec<u32> = times.iter().map(|&x| (x * 1000.0) as u32).collect();
    Ok(Arc::from(times))
}

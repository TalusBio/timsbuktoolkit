use clap::Parser;
use index::InputQuery;
use serde_json::{
    Value,
    json,
};
use std::io::{
    BufReader,
    Read,
    Write,
};
use std::net::{
    TcpListener,
    TcpStream,
};
use std::sync::Arc;
use std::thread;
use timsseek::errors::{
    Result,
    TimsSeekError,
};
use timsseek::scoring::full_results::FullQueryResult;

mod cli;
mod index;

struct DaemonServer {
    index: Arc<index::BundledDotDIndex>,
    running: std::sync::atomic::AtomicBool,
}

impl DaemonServer {
    pub fn new(index: index::BundledDotDIndex) -> std::io::Result<Self> {
        Ok(Self {
            index: Arc::new(index),
            running: std::sync::atomic::AtomicBool::new(true),
        })
    }

    pub fn run(&self, addr: &str) -> std::io::Result<()> {
        let listener = TcpListener::bind(addr)?;
        println!("Listening on {}", addr);

        while self.running.load(std::sync::atomic::Ordering::Relaxed) {
            // listener.set_nonblocking(true)?;

            match listener.accept() {
                Ok((stream, _)) => {
                    let index = Arc::clone(&self.index);
                    let running = &self.running;

                    match handle_connection(stream, index, running) {
                        Ok(_) => (),
                        Err(e) => eprintln!("Error handling connection: {}", e),
                    };
                }
                Err(ref e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                    thread::sleep(std::time::Duration::from_millis(100));
                    continue;
                }
                Err(e) => eprintln!("Error accepting connection: {}", e),
            }
        }

        Ok(())
    }
}

fn handle_connection(
    mut stream: TcpStream,
    index: Arc<index::BundledDotDIndex>,
    _running: &std::sync::atomic::AtomicBool,
) -> std::io::Result<()> {
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut buffer = String::new();

    loop {
        buffer.clear();
        reader.read_to_string(&mut buffer)?;

        if buffer.is_empty() {
            break;
        }
        println!("read data: {}", buffer);

        let query: Value = match serde_json::from_str(&buffer) {
            Ok(q) => q,
            Err(e) => {
                let response = json!({
                    "status": "error",
                    "data": format!("Invalid JSON format: {}", e)
                });
                send_response(&mut stream, &response)?;
                continue;
            }
        };

        let query: InputQuery = match serde_json::from_value(query) {
            Ok(q) => q,
            Err(e) => {
                let response = json!({
                    "status": "error",
                    "data": format!("Invalid query format: {}", e)
                });
                send_response(&mut stream, &response)?;
                continue;
            }
        };

        let start = std::time::Instant::now();
        let query_res: Result<FullQueryResult> = index.query(query.into());
        let elap_time = start.elapsed();
        println!("Querying took {:#?} for query", elap_time);
        let response = match query_res {
            Ok(q) => json!({
                "status": "success",
                "data": q
            }),
            Err(e) => json!({
                "status": "error",
                "data": format!("{}", e)
            }),
        };
        send_response(&mut stream, &response)?;
    }

    Ok(())
}

fn send_response(stream: &mut TcpStream, response: &Value) -> std::io::Result<()> {
    stream.write_all(response.to_string().as_bytes())?;
    stream.write_all(b"\n")?;
    Ok(())
}

// Example usage
fn main() -> Result<()> {
    let conf = cli::Cli::parse();
    let sample = InputQuery::sample();
    let tol = conf.read_config()?;
    let index = index::BundledDotDIndex::new(conf.dotd_file, tol)?;

    println!("Starting server");
    println!(
        "Sample query: \n{}",
        serde_json::to_string_pretty(&sample).unwrap()
    );

    let st = std::time::Instant::now();
    let _check_out = match index.query(sample.into()) {
        Ok(q) => {
            println!("Query OK");
            q
        }
        Err(e) => {
            println!("Query failed: {}", e);
            return Err(e);
        }
    };
    let elap_time = st.elapsed();
    println!("Querying took {:#?} for sample query", elap_time);

    // println!("Query result: \n{}", serde_json::to_string_pretty(&check_out).unwrap());

    let server = match DaemonServer::new(index) {
        Ok(s) => s,
        Err(e) => {
            return Err(TimsSeekError::Io {
                source: e,
                path: None,
            });
        }
    };
    match server.run(&conf.address) {
        Ok(_) => Ok(()),
        Err(e) => Err(TimsSeekError::Io {
            source: e,
            path: None,
        }),
    }
}

// Example client code
// fn example_client() -> Result<()> {
//     let mut stream = TcpStream::connect("127.0.0.1:8080")?;
//
//     let query = json!({
//         "type": "search",
//         "term": "example"
//     });
//
//     match stream.write_all(query.to_string().as_bytes()) {
//         Ok(_) => (),
//         Err(e) => {
//             return Err(TimsSeekError::Io {
//                 source: e,
//                 path: None,
//             });
//         }
//     };
//     match stream.write_all(b"\n") {
//         Ok(_) => (),
//         Err(e) => {
//             return Err(TimsSeekError::Io {
//                 source: e,
//                 path: None,
//             });
//         }
//     };
//
//     let mut reader = BufReader::new(stream);
//     let mut response = String::new();
//     reader.read_line(&mut response)?;
//
//     let response_json: Value = serde_json::from_str(&response).unwrap();
//     println!("Response: {}", response_json);
//
//     Ok(())
// }

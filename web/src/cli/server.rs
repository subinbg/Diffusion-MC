//! Development HTTP server for serving the site.

use std::fs;
use std::io::Read;
use std::path::Path;
use tiny_http::{Response, Server};

/// Serve the site with a simple HTTP server.
pub fn serve_site(port: u16) -> Result<(), Box<dyn std::error::Error>> {
    let addr = format!("0.0.0.0:{}", port);
    let server = Server::http(&addr).map_err(|e| format!("Failed to start server: {}", e))?;

    println!("Serving site at http://localhost:{}", port);
    println!("Press Ctrl+C to stop");

    for request in server.incoming_requests() {
        let url = request.url().to_string();
        let path = if url == "/" {
            "dist/index.html".to_string()
        } else if url.ends_with('/') {
            format!("dist{}index.html", url)
        } else {
            format!("dist{}", url)
        };

        let response = match serve_file(&path) {
            Ok(resp) => resp,
            Err(_) => {
                // Try with .html extension
                match serve_file(&format!("{}.html", path)) {
                    Ok(resp) => resp,
                    Err(_) => {
                        let not_found = "404 Not Found";
                        Response::from_string(not_found).with_status_code(404)
                    }
                }
            }
        };

        let _ = request.respond(response);
    }

    Ok(())
}

/// Serve a single file.
fn serve_file(path: &str) -> Result<Response<std::io::Cursor<Vec<u8>>>, std::io::Error> {
    let mut file = fs::File::open(path)?;
    let mut content = Vec::new();
    file.read_to_end(&mut content)?;

    let content_type = match Path::new(path).extension().and_then(|e| e.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        Some("js") => "application/javascript; charset=utf-8",
        Some("wasm") => "application/wasm",
        Some("json") => "application/json; charset=utf-8",
        Some("png") => "image/png",
        Some("jpg") | Some("jpeg") => "image/jpeg",
        Some("svg") => "image/svg+xml",
        _ => "application/octet-stream",
    };

    Ok(Response::from_data(content)
        .with_header(tiny_http::Header::from_bytes("Content-Type", content_type).unwrap())
        .with_header(
            tiny_http::Header::from_bytes("Cross-Origin-Opener-Policy", "same-origin").unwrap(),
        )
        .with_header(
            tiny_http::Header::from_bytes("Cross-Origin-Embedder-Policy", "require-corp").unwrap(),
        ))
}

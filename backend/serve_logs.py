from fastapi import FastAPI
from fastapi.responses import HTMLResponse
import os
import uvicorn

app = FastAPI()

LOG_FILE = "ambeed_crawl.log"

@app.get("/", response_class=HTMLResponse)
def read_logs():
    if not os.path.exists(LOG_FILE):
        return "Log file not found yet. Crawler might be initializing..."
    
    try:
        # Read last 300 lines safely
        with open(LOG_FILE, "r", encoding="utf-8", errors='replace') as f:
            lines = f.readlines()
            last_lines = lines[-300:] if len(lines) > 300 else lines
        
        content = "".join(last_lines)
    except Exception as e:
        content = f"Error reading log: {e}"
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
        <head>
            <title>🕷️ Ambeed Crawler Logs</title>
            <meta http-equiv="refresh" content="2"> <!-- Auto Refresh every 2s -->
            <style>
                body {{ 
                    background-color: #1e1e1e; 
                    color: #00ff00; 
                    font-family: 'Consolas', 'Courier New', monospace; 
                    padding: 20px;
                    font-size: 14px;
                }}
                h1 {{ border-bottom: 1px solid #444; padding-bottom: 10px; }}
                pre {{ white-space: pre-wrap; word-wrap: break-word; }}
                .status {{ color: #ffff00; font-weight: bold; }}
            </style>
        </head>
        <body>
            <h1>🚀 Ambeed Crawler Real-time Logs</h1>
            <div class="status">Monitoring: {LOG_FILE} (Auto-refresh: 2s)</div>
            <pre>{content}</pre>
            <script>
                window.scrollTo(0, document.body.scrollHeight);
            </script>
        </body>
    </html>
    """
    return html_content

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8001)

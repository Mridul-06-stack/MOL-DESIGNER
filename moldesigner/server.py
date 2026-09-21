"""Server module – FastAPI + WebSocket / SSE entry point."""

import uvicorn


def run_server(host: str = "127.0.0.1", port: int = 8000, reload: bool = False):
    """Run the MolDesigner API server."""
    uvicorn.run("server.api:app", host=host, port=port, reload=reload)


if __name__ == "__main__":
    run_server()

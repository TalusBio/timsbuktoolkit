import json
import socket

from .constants import DEFAULT_PORT


def query_server(host="localhost", port: int = DEFAULT_PORT, query_data: dict = {}):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.connect((host, port))

        # Send query
        query = json.dumps(query_data).encode("utf-8")
        s.sendall(query)
        s.shutdown(socket.SHUT_WR)

        # Read response with timeout
        s.settimeout(2.0)  # 2 second timeout
        chunks = []

        try:
            while True:
                chunk = s.recv(2048)
                if not chunk:
                    break
                chunks.append(chunk)
        except socket.timeout:
            pass  # It's okay if we timeout after receiving data

        # print(chunks)
        response = b"".join(chunks).decode("utf-8")

        try:
            return json.loads(response)
        except json.JSONDecodeError:
            return response

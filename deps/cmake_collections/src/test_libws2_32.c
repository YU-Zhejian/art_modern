#include <winsock2.h>
#include <ws2tcpip.h>
#include <stdio.h>

int main() {
    WSADATA wsaData;
    int result = WSAStartup(MAKEWORD(2,2), &wsaData);
    if (result != 0) {
        return 1;
    }
    SOCKET sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (sock == INVALID_SOCKET) {
        WSACleanup();
        return 1;
    }
    closesocket(sock);
    WSACleanup();
    return 0;
}

#pragma once

class FileStreamInterface {
public:
    virtual void write(const std::string& str) { }

    virtual void close() { }
};

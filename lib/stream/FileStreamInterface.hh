#pragma once

class FileStreamInterface {
public:
    virtual void write(const std::string& str)
    {
        // Implementation-defined
    }

    virtual void close()
    {
        // Implementation-defined
    }

    virtual ~FileStreamInterface() = default;
};

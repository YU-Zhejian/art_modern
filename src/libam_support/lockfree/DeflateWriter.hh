#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <libdeflate.h>

#include <zlib.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class WriterInterface {
public:
    WriterInterface() = default;
    virtual ~WriterInterface() = default;
    virtual void write(std::string&& data) = 0;
    virtual void write(const std::string& data) = 0;
    DELETE_COPY(WriterInterface)
    DELETE_MOVE(WriterInterface)
    virtual void close() = 0;
    virtual void flush() = 0;
    [[nodiscard]] virtual std::size_t get_total_bytes_written() const = 0;
};

class SimpleWriter : public WriterInterface {
public:
    DELETE_COPY(SimpleWriter)
    DELETE_MOVE(SimpleWriter)
    explicit SimpleWriter(const std::string& filename)
        : ofs_(filename, std::ios::binary)
    {
    }
    ~SimpleWriter() override { SimpleWriter::close(); }

    void write(std::string&& data) override
    {
        std::string const data_in = std::move(data);
        ofs_ << data_in;
        total_bytes_written_ += data_in.size();
    }
    void close() override
    {
        if (ofs_.is_open()) {
            ofs_.close();
        }
    }
    void flush() override { ofs_.flush(); }
    std::size_t get_total_bytes_written() const override { return total_bytes_written_; }
    void write(const std::string& data) override
    {
        ofs_ << data;
        total_bytes_written_ += data.size();
    }

private:
    std::ofstream ofs_;
    std::size_t total_bytes_written_ { 0 };
};

class ZLibWriter final : public WriterInterface {
public:
    explicit ZLibWriter(const std::string& filename, const int level = 6, const std::size_t buffer_size = 1 << 20)
        : buffer_size_(buffer_size)
    {
        gz_file_ = gzopen(filename.c_str(), ("wb" + std::to_string(level)).c_str());
        if (gz_file_ == nullptr) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        gzbuffer(gz_file_, static_cast<unsigned int>(buffer_size_));
    }

    ~ZLibWriter() override
    {
        ZLibWriter::close();
        gzclose(gz_file_);
    }
    DELETE_COPY(ZLibWriter)
    DELETE_MOVE(ZLibWriter)

    void write(std::string&& data) override
    {
        std::string const data_in = std::move(data);
        auto bytes_written = gzwrite(gz_file_, data_in.data(), static_cast<unsigned int>(data_in.size()));
        if (bytes_written < 0) {
            throw std::runtime_error("Failed to write to gzip file.");
        }
        update_offset_();
    }
    void write(const std::string& data) override
    {
        auto bytes_written = gzwrite(gz_file_, data.data(), static_cast<unsigned int>(data.size()));
        if (bytes_written < 0) {
            throw std::runtime_error("Failed to write to gzip file.");
        }
        update_offset_();
    }

    void flush() override
    {
        gzflush(gz_file_, Z_SYNC_FLUSH);
        update_offset_();
    }

    void close() override
    {
        if (gz_file_ != nullptr) {
            update_offset_();
            gzclose(gz_file_);
            gz_file_ = nullptr;
        }
    }

    [[nodiscard]] std::size_t get_total_bytes_written() const override { return total_bytes_written_; }

private:
    gzFile gz_file_;
    std::size_t buffer_size_;
    std::size_t total_bytes_written_ { 0 };
    void update_offset_()
    {
        auto current_offset =
#ifdef Z_LARGE64
            gzoffset64
#else
            gzoffset
#endif
            (gz_file_);
        if (current_offset < 0) {
            throw std::runtime_error("Failed to get current offset in gzip file.");
        }
        total_bytes_written_ = current_offset;
    }
};

class DeflateWriter final : public WriterInterface {
public:
    // Compression levels usually range from 1 to 12
    explicit DeflateWriter(const std::string& filename, const int level = 6, const std::size_t buffer_size = 1 << 20)
        : compressor_(libdeflate_alloc_compressor(level))
        , buffer_size_(buffer_size)
    {
        out_stream_.open(filename, std::ios::out | std::ios::binary);
        buffer_.reserve(buffer_size_);
    }

    ~DeflateWriter() override
    {
        DeflateWriter::close();
        libdeflate_free_compressor(compressor_);
    }

    DELETE_COPY(DeflateWriter)
    DELETE_MOVE(DeflateWriter)

    // Write moves the string into our internal buffer
    void write(std::string&& data) override
    {
        const std::string data_in = std::move(data);
        buffer_.insert(buffer_.end(), data_in.begin(), data_in.end());
        if (buffer_.size() >= buffer_size_) {
            flush_buffer_to_disk();
        }
    }
    void write(const std::string& data) override
    {
        buffer_.insert(buffer_.end(), data.begin(), data.end());
        if (buffer_.size() >= buffer_size_) {
            flush_buffer_to_disk();
        }
    }

    // Manual flush: Compresses whatever is currently in the buffer
    void flush() override
    {
        if (!buffer_.empty()) {
            flush_buffer_to_disk();
        }
        out_stream_.flush();
    }

    void close() override
    {
        if (out_stream_.is_open()) {
            flush();
            out_stream_.close();
        }
    }

    std::size_t get_total_bytes_written() const override { return compressed_size_; }

private:
    void flush_buffer_to_disk()
    {
        if (buffer_.empty()) {
            return;
        }

        // libdeflate needs a destination buffer.
        // Worst case: compressed size slightly larger than input
        const size_t max_compressed_size = libdeflate_gzip_compress_bound(compressor_, buffer_.size());
        compressed_out_.resize(max_compressed_size);

        const size_t actual_size = libdeflate_gzip_compress(
            compressor_, buffer_.data(), buffer_.size(), compressed_out_.data(), compressed_out_.size());

        if (actual_size > 0) {
            // Write the compressed size first (common for block-based custom formats)
            // or just write the raw block if you're using a specific container format.
            out_stream_.write(compressed_out_.data(), static_cast<std::streamsize>(actual_size));
        }
        compressed_size_ += actual_size;
        buffer_.clear();
    }

    std::ofstream out_stream_;
    std::vector<char> buffer_;
    libdeflate_compressor* compressor_;
    std::size_t compressed_size_ = 0;
    std::vector<char> compressed_out_;
    std::size_t buffer_size_;
};

#pragma once

#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace labw::art_modern {
class GTfAttributesTokenizerException : public std::runtime_error {
public:
    explicit GTfAttributesTokenizerException(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

class GtfAttributesTokens {
public:
    std::vector<std::string> keys;
    std::vector<std::string> values;
};

class GtfAttributesTokenizerWithFSA {
public:
    GtfAttributesTokenizerWithFSA();
    GtfAttributesTokens parse(const std::string& attributes_to_parse);

private:
    enum class LexerStatus : std::uint8_t { WAITING_FOR_KEY, EXTENDING_KEY, WAITING_FOR_VALUE, EXTENDING_VALUE };
    enum class QuotationMarkStatus : std::uint8_t { UNSET, SET_LEFT_SINGLE, SET_LEFT_DOUBLE };
    enum class TokenType : std::uint8_t { KEY, VALUE };
    std::vector<std::string> tokens_for_keys_;
    std::vector<std::string> tokens_for_values_;
    std::string input_attributes_str_;
    std::size_t current_lexer_position_ = 0;
    std::string current_token_value_;
    LexerStatus current_lexer_status_ = LexerStatus::WAITING_FOR_KEY;
    TokenType current_token_type_ = TokenType::VALUE;
    TokenType last_token_type_ = TokenType::VALUE;
    std::size_t quotation_start_ = 0;
    QuotationMarkStatus quotation_mark_status_ = QuotationMarkStatus::UNSET;

    void reset_();
    void reset_quotation_mark_status_();
    static QuotationMarkStatus is_quote_(char c);
    [[nodiscard]] bool is_peekable_() const;
    [[nodiscard]] char peek_() const;
    [[nodiscard]] char get_() const;
    void add_parsed_token_to_token_list_();
    void extend_current_token_(char getc);
    void throw_waiting_for_value_() const;
    void throw_waiting_for_end_quote_() const;
    void case_wait_for_key_(char getc);
    void case_wait_for_value_(char getc);
    void case_extending_key_(char getc);
    void case_extending_value_(char getc);
};

} // namespace labw::art_modern

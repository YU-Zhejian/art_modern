#include "libam_support/ref/gtf/GtfLex.hh"

#include <cstdlib>
#include <string>
#include <vector>

namespace labw::art_modern {
namespace {
    inline bool isResp_(const char c) { return c == ';'; }
    inline bool isBlank_(const char c) { return c == ' ' || c == '\t'; }
} // namespace

GtfAttributesTokenizerWithFSA::GtfAttributesTokenizerWithFSA() { reset_(); }

void GtfAttributesTokenizerWithFSA::reset_()
{
    tokens_for_keys_.clear();
    tokens_for_values_.clear();
    current_lexer_position_ = 0;
    current_token_value_.clear();
    current_lexer_status_ = LexerStatus::WAITING_FOR_KEY;
    last_token_type_ = TokenType::VALUE;
    reset_quotation_mark_status_();
}

void GtfAttributesTokenizerWithFSA::reset_quotation_mark_status_()
{
    quotation_start_ = 0;
    quotation_mark_status_ = QuotationMarkStatus::UNSET;
}

GtfAttributesTokenizerWithFSA::QuotationMarkStatus GtfAttributesTokenizerWithFSA::is_quote_(const char c)
{
    switch (c) {
    case '\'':
        return QuotationMarkStatus::SET_LEFT_SINGLE;
    case '\"':
        return QuotationMarkStatus::SET_LEFT_DOUBLE;
    default:
        return QuotationMarkStatus::UNSET;
    }
}

bool GtfAttributesTokenizerWithFSA::is_peekable_() const
{
    return current_lexer_position_ + 1 < input_attributes_str_.length();
}

char GtfAttributesTokenizerWithFSA::peek_() const { return input_attributes_str_.at(current_lexer_position_ + 1); }

char GtfAttributesTokenizerWithFSA::get_() const { return input_attributes_str_.at(current_lexer_position_); }

void GtfAttributesTokenizerWithFSA::add_parsed_token_to_token_list_()
{
    if (last_token_type_ == current_token_type_) {
        throw GTfAttributesTokenizerException("Get two tokens of a kind at " + std::to_string(current_lexer_position_));
    }
    if (current_token_type_ == TokenType::KEY) {
        tokens_for_keys_.push_back(current_token_value_);
    } else {
        tokens_for_values_.push_back(current_token_value_);
    }
    last_token_type_ = current_token_type_;
}

void GtfAttributesTokenizerWithFSA::extend_current_token_(const char getc) { current_token_value_.push_back(getc); }

void GtfAttributesTokenizerWithFSA::throw_waiting_for_value_() const
{
    throw GTfAttributesTokenizerException("Waiting for values at " + std::to_string(current_lexer_position_));
}

void GtfAttributesTokenizerWithFSA::throw_waiting_for_end_quote_() const
{
    throw GTfAttributesTokenizerException("Waiting for end quote ("
        + std::to_string(static_cast<int>(quotation_mark_status_)) + ") at " + std::to_string(current_lexer_position_)
        + " started at " + std::to_string(quotation_start_));
}

void GtfAttributesTokenizerWithFSA::case_wait_for_key_(const char getc)
{
    if (!(isBlank_(getc) || isResp_(getc))) {
        if (is_quote_(getc) != QuotationMarkStatus::UNSET) {
            throw GTfAttributesTokenizerException(
                "Error at " + std::to_string(current_lexer_position_) + ": Quotes not allowd in keys.");
        }
        current_lexer_status_ = LexerStatus::EXTENDING_KEY;
        current_token_value_.clear();
        current_token_type_ = TokenType::KEY;
        extend_current_token_(getc);
    }
}

void GtfAttributesTokenizerWithFSA::case_wait_for_value_(const char getc)
{
    if (!isBlank_(getc)) {
        if (is_quote_(getc) != QuotationMarkStatus::UNSET) {
            quotation_mark_status_ = is_quote_(getc);
            quotation_start_ = current_lexer_position_;
            current_lexer_status_ = LexerStatus::EXTENDING_VALUE;
            current_token_value_.clear();
            current_token_type_ = TokenType::VALUE;
        } else {
            current_lexer_status_ = LexerStatus::EXTENDING_VALUE;
            current_token_value_.clear();
            current_token_type_ = TokenType::VALUE;
            extend_current_token_(getc);
        }
    }
}

void GtfAttributesTokenizerWithFSA::case_extending_key_(const char getc)
{
    if (isBlank_(getc)) {
        current_lexer_status_ = LexerStatus::WAITING_FOR_VALUE;
        add_parsed_token_to_token_list_();
    } else {
        extend_current_token_(getc);
    }
}

void GtfAttributesTokenizerWithFSA::case_extending_value_(const char getc)
{
    if (quotation_mark_status_ != QuotationMarkStatus::UNSET) {
        if (quotation_mark_status_ == is_quote_(getc)) {
            if (is_peekable_()) {
                char const peeked = peek_();
                if (isBlank_(peeked) || isResp_(peeked)) {
                    current_lexer_status_ = LexerStatus::WAITING_FOR_KEY;
                    reset_quotation_mark_status_();
                    add_parsed_token_to_token_list_();
                } else {
                    throw GTfAttributesTokenizerException("Unexpected non-whitespace/termination character ("
                        + std::to_string(peeked) + ") after termination of quote at "
                        + std::to_string(current_lexer_position_));
                }
            } else {
                current_lexer_status_ = LexerStatus::WAITING_FOR_KEY;
                reset_quotation_mark_status_();
                add_parsed_token_to_token_list_();
            }
        } else {
            extend_current_token_(getc);
        }
    } else if (isBlank_(getc) || isResp_(getc)) {
        current_lexer_status_ = LexerStatus::WAITING_FOR_KEY;
        reset_quotation_mark_status_();
        add_parsed_token_to_token_list_();
    } else {
        extend_current_token_(getc);
    }
}

GtfAttributesTokens GtfAttributesTokenizerWithFSA::parse(const std::string& attributes_to_parse)
{
    reset_();
    this->input_attributes_str_ = attributes_to_parse;

    char getc = 0;
    while (current_lexer_position_ < input_attributes_str_.length()) {
        getc = get_();
        switch (current_lexer_status_) {
        case LexerStatus::WAITING_FOR_KEY:
            case_wait_for_key_(getc);
            break;
        case LexerStatus::EXTENDING_KEY:
            case_extending_key_(getc);
            break;
        case LexerStatus::WAITING_FOR_VALUE:
            case_wait_for_value_(getc);
            break;
        case LexerStatus::EXTENDING_VALUE:
            case_extending_value_(getc);
            break;
        }
        current_lexer_position_ += 1;
    }
    switch (current_lexer_status_) {
    case LexerStatus::WAITING_FOR_KEY:
        // The normal situation
        break;
    case LexerStatus::EXTENDING_KEY:
        if (quotation_mark_status_ == QuotationMarkStatus::UNSET) {
            throw_waiting_for_value_();
        } else {
            throw_waiting_for_end_quote_();
        }
        break;
    case LexerStatus::WAITING_FOR_VALUE:
        throw_waiting_for_value_();
        break;
    case LexerStatus::EXTENDING_VALUE:
        if (quotation_mark_status_ == QuotationMarkStatus::UNSET) {
            add_parsed_token_to_token_list_();
        } else {
            throw_waiting_for_end_quote_();
        }
        break;
    }
    if (tokens_for_keys_.size() != tokens_for_values_.size()) {
        throw GTfAttributesTokenizerException("Different number of keys (" + std::to_string(tokens_for_keys_.size())
            + ") and values (" + std::to_string(tokens_for_values_.size()) + ") at "
            + std::to_string(current_lexer_position_));
    }
    GtfAttributesTokens result;
    result.keys = tokens_for_keys_;
    result.values = tokens_for_values_;
    return result;
}
} // namespace labw::art_modern

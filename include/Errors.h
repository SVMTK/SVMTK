#ifndef  __Errors_H

#define __Errors_H



class PreconditionError: public std::exception {
public:
    explicit PreconditionError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
} ;

class EmptyMeshError: public std::exception {
public:
    explicit EmptyMeshError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
} ;


class InvalidArgumentError: public std::exception {
public:
    explicit InvalidArgumentError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
} ;


class AlgorithmError: public std::exception {
public:
    explicit AlgorithmError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
} ;


#endif

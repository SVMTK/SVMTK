#ifndef  __Errors_H

#define __Errors_H
#include<string>
/**
 * \class PreconditionError
 *
 * Handles the case of precondition error in SVMTK.
 *
 * @throws std::exception 
 */ 
class PreconditionError: public std::exception {
public:
    explicit PreconditionError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};

/**
 * \class EmptyMeshError
 *
 * Handles the case of a empty mesh in SVMTK class object.
 *
 * @throws std::exception 
 */ 
class EmptyMeshError: public std::exception {
public:
    explicit EmptyMeshError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};

/**
 * \class InvalidArgumentError
 *
 * Handles the case of invalidArgument error in SVMTK.
 *
 * @throws std::exception 
 */ 
class InvalidArgumentError: public std::exception {
public:
    explicit InvalidArgumentError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};

/**
 * \class AlgorithmError
 *
 * Handles the case of algorithm error in SVMTK.
 *
 * @throws std::exception 
 */ 
class AlgorithmError: public std::exception {
public:
    explicit AlgorithmError(const char * m) : message{m} {}
    const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};


#endif

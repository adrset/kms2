#include "Properties.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
Properties::Properties(std::string name) : m_fileName (name){

    std::string attributeName;
    double attributeValue;
    std::ifstream propFile(m_fileName);

    while(propFile>>attributeName>>attributeValue){

        std::cout<<"Attribute:" << attributeName << ", value: " << attributeValue << std::endl;
        if (m_values.find(attributeName) == m_values.end()) {
            m_values.insert(std::make_pair(attributeName , attributeValue ));
        }
        else {
           throw std::runtime_error("Error!\nAttribute '" + attributeName + "' was already defined!"); 
        }
    }

}

double Properties::getProperty(std::string prop){
    if (m_values.find(prop) != m_values.end()) {
        return m_values[prop];
    }else{
        throw std::runtime_error("Error!\nAttribute '" + prop + "' not defined!"); 
    }

}
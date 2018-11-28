#ifndef PROPERTIES_H_
#define PROPERTIES_H_
#include <iostream>
#include <map>
	class Properties{
		public:

			Properties(std::string);
			double getProperty(std::string);

		private:
			std::string m_fileName;
			std::map<std::string, double> m_values;

	};

#endif

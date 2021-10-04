#ifndef CHE_OFF_H
#define CHE_OFF_H

#include "che.h"


class che_off : public che {
    public:
        explicit che_off(const std::string &file);
        virtual ~che_off();
    private:
        void read_file(const std::string &file);
};




#endif // CHE_OFF_H


#include "Json/json.hpp"
#include "config.hpp"

class Context : public config_t
{
    private:
        nlohmann::json data_;
        void merge_dict(data_, input_schema);

}





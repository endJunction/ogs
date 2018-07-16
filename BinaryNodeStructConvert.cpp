#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <iterator>

struct Node
{
    long id;
    double x;
    double y;
    double z;
};

int main(int argc, char* argv[])
{
    std::string const filename = argv[2];
    std::ifstream in(filename, std::ios::binary);
    if (!in)
        return EXIT_FAILURE;
    while (!in.eof())
    {
        if (std::string(argv[1]) == "n")
        {
            Node node;
            in.read(reinterpret_cast<char*>(&node), sizeof(Node));
            std::cout << node.id << " : (" << node.x << ", " << node.y << ", "
                      << node.z << "); ";
        }
        else if (std::string(argv[1]) == "c")
        {
            long value[14];
            if (in.read(reinterpret_cast<char*>(&value), 14 * sizeof(long)))
            {
                std::cout << "nodes.size() : " << value[0] << "\n"
                          << "number_of_base_nodes), : " << value[1] << "\n"
                          << "regular_elements.size()), : " << value[2] << "\n"
                          << "ghost_elements.size()), : " << value[3] << "\n"
                          << "number_of_non_ghost_base_nodes), : " << value[4]
                          << "\n"
                          << "number_of_non_ghost_nodes), : " << value[5]
                          << "\n"
                          << "number_of_mesh_base_nodes), : " << value[6]
                          << "\n"
                          << "number_of_mesh_all_nodes), : " << value[7] << "\n"
                          << "getNumberOfIntegerVariablesOfElements(regular_"
                             "elements)), "
                             ": "
                          << value[8] << "\n"
                          << "getNumberOfIntegerVariablesOfElements(ghost_"
                             "elements)), : "
                          << value[9] << "\n"
                          << "offsets(0), : " << value[10] << "\n"
                          << "offsets(1), : " << value[11] << "\n"
                          << "offsets(2), : " << value[12] << "\n"
                          << "reserved, : " << value[13] << "\n\n";
            }
            else
            {
                std::cerr << "Could not read 14*long";
            }
        }
        else
        {
            long value;
            in.read(reinterpret_cast<char*>(&value), sizeof(long));
            std::cout << value << " ";
        }
    }
    return EXIT_SUCCESS;
}

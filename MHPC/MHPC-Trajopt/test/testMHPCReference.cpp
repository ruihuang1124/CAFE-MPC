#include "MHPCReference.h"
#include "QuadReference.h"

int main()
{
    WBReference wb_ref;
    SRBReference srb_ref;

    std::string reference_file_path = "../Reference/Data/quad_reference.csv";       
    std::shared_ptr<QuadReference> quad_ref;
    quad_ref = std::make_shared<QuadReference>();
    quad_ref->load_top_level_data(reference_file_path);

    wb_ref.set_quadruped_reference(quad_ref);
    srb_ref.set_quadruped_reference(quad_ref);
}
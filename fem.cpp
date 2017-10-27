/*! @file

    @brief A quick look into finite elements
*/

#include <stdio.h>
#include <assert.h>

#include "linalgcpp.hpp"

using namespace linalgcpp;

SparseMatrix<double> ElemDataDirect()
{
    SparseMatrix<int> elem_node = ReadTable("element_node_4.txt");
    std::vector<double> node_coo = ReadText<double>("coordinati_4.txt");
    std::vector<double> elem_vect = ReadText<double>("element_matrix_4.txt");

    const size_t elem_size = 3;
    const size_t num_elems = elem_node.Rows();
    const size_t num_nodes = elem_node.Cols();

    assert(node_coo.size() == 2 * num_nodes);
    assert(elem_vect.size() == elem_size * elem_size * num_elems);

    CooMatrix<double> coo_assembly;

    size_t count = 0;

    for (size_t elem = 0; elem < num_elems; ++elem)
    {
        const auto indices = elem_node.GetIndices(elem);

        for (size_t i = 0; i < elem_size; ++i)
        {
            for (size_t j = 0; j < elem_size; ++j)
            {
                const size_t x = indices[i];
                const size_t y = indices[j];
                const double val = elem_vect[count++];

                coo_assembly.Add(x, y, val);
            }
        }
    }

    return coo_assembly.ToSparse();
}

template <typename T = int>
SparseMatrix<T> MakeEdgeVertex(const SparseMatrix<T>& mat)
{
    assert(mat.Rows() == mat.Cols());
    assert((mat.nnz() - mat.Rows()) % 2 == 0);

    const size_t num_vertices = mat.Rows();
    const size_t num_edges = (mat.nnz() - mat.Rows()) / 2;

    std::vector<int> indptr(num_edges + 1);
    std::vector<int> indices(num_edges * 2);
    std::vector<T> data(num_edges * 2, 1);

    indptr[0] = 0;

    int count = 0;

    for (size_t i = 0; i < num_vertices; ++i)
    {
        for (int j = mat.GetIndptr()[i]; j < mat.GetIndptr()[i + 1]; ++j)
        {
            const size_t col = mat.GetIndices()[j];

            if (i < col)
            {
                indices[count++] = i;
                indices[count++] = col;

                indptr[count / 2] = count;
            }
        }
    }

    return SparseMatrix<T>(indptr, indices, data, num_edges, num_vertices);
}

template <typename T = int>
SparseMatrix<T> RestrictInterior(const SparseMatrix<T>& mat)
{
    std::vector<int> indptr(mat.Rows() + 1);
    std::vector<int> indices;

    indptr[0] = 0;

    for (size_t i = 0; i < mat.Rows(); ++i)
    {
        for (int j = mat.GetIndptr()[i]; j < mat.GetIndptr()[i + 1]; ++j)
        {
            if (mat.GetData()[j] >= 2)
            {
                indices.push_back(mat.GetIndices()[j]);
            }
        }

        indptr[i + 1] = indices.size();
    }

    std::vector<T> data(indices.size(), 1);

    return SparseMatrix<T>(indptr, indices, data, mat.Rows(), mat.Cols());
}

SparseMatrix<double> ElemData()
{
    std::vector<double> elem_vect = ReadText<double>("element_matrix_4.txt");
    std::vector<double> node_coo = ReadText<double>("coordinati_4.txt");
    SparseMatrix<int> elem_node = ReadTable("element_node_4.txt");

    const size_t elem_size = 3 * 3;
    const size_t num_elems = elem_node.Rows();
    const size_t num_nodes = elem_node.Cols();

    assert(node_coo.size() == 2 * num_nodes);
    assert(elem_vect.size() == elem_size * num_elems);

    CooMatrix<double> coo_assembly;

    auto data_ptr = std::begin(elem_vect);

    for (size_t elem = 0; elem < num_elems; ++elem)
    {
        printf("Element: %ld\nNodes:\n", elem);

        const auto indices = elem_node.GetIndices(elem);

        for (auto index : indices)
        {
            const double x_i = node_coo[index * 2];
            const double y_i = node_coo[index * 2 + 1];

            printf("\t%d:\t(%.3f, %.3f)\n", index, x_i, y_i);
        }

        std::vector<double> elem_data(data_ptr, data_ptr + elem_size);
        data_ptr += elem_size;

        DenseMatrix element(3, 3, elem_data);

        element.Print("Elem Mat:");

        coo_assembly.Add(indices, element);

        printf("------------------------------\n\n");
    }

    //coo_assembly.Print("coo:");

    return coo_assembly.ToSparse();
}

void FindBoundary(const SparseMatrix<int>& elem_node,
                  const std::vector<double>& node_coo)
{
    // Create relationships
    SparseMatrix<int> node_elem = elem_node.Transpose();
    SparseMatrix<int> node_node = node_elem.Mult(elem_node);

    SparseMatrix<int> edge_node = MakeEdgeVertex(node_node);
    SparseMatrix<int> node_edge = edge_node.Transpose();

    SparseMatrix<int> elem_edge = RestrictInterior(elem_node.Mult(node_edge));

    SparseMatrix<int> edge_elem = elem_edge.Transpose();

    // Mark boundaries
    const size_t num_elems = elem_node.Rows();
    const size_t num_nodes = node_elem.Rows();
    const size_t num_edges = edge_elem.Rows();

    std::vector<int> elem_boundary(num_elems, 0);
    std::vector<int> node_boundary(num_nodes, 0);
    std::vector<int> edge_boundary(num_edges, 0);

    for (size_t edge = 0; edge < num_edges; ++edge)
    {
        assert(edge_elem.RowSize(edge) > 0);
        assert(edge_elem.RowSize(edge) <= 2);

        if (edge_elem.RowSize(edge) == 1)
        {
            edge_boundary[edge]++;

            for (auto node : edge_node.GetIndices(edge))
            {
                node_boundary[node]++;
            }

            for (auto elem : edge_elem.GetIndices(edge))
            {
                elem_boundary[elem]++;
            }
        }
    }

    // Check if all boundaries are correct
    for (size_t node = 0; node < num_nodes; ++node)
    {
        double x_i = node_coo[node * 2];
        double y_i = node_coo[node * 2 + 1];

        if (node_boundary[node] > 0)
        {
            //printf("%ld (%.2f, %.2f)\n", node, x_i, y_i);
            assert(x_i == 0 || x_i == 1 || y_i == 0 || y_i == 1);
        }
        else
        {
            assert(x_i != 0 && x_i != 1 && y_i != 0 && y_i != 1);
        }
    }
}

int main(int argc, char** argv)
{
    auto sp_elem = ElemData();
    auto sp_direct = ElemDataDirect();

    // Read in some data
    SparseMatrix<int> elem_node = ReadTable("element_node_4.txt");
    std::vector<double> node_coo = ReadText<double>("coordinati_4.txt");

    FindBoundary(elem_node, node_coo);

    assert(sp_elem.Rows() == sp_elem.Cols());
    assert(sp_direct.Rows() == sp_direct.Cols());
    assert(sp_elem.Rows() == sp_direct.Rows());

    auto size = sp_elem.Rows();

    auto x = Vector<double>(size);
    auto y = Vector<double>(size);

    Randomize(y, -1, 1);
    Randomize(x, -1, 1);

    auto sp_x = sp_elem.Mult(x);
    auto sp_direct_x = sp_direct.Mult(x);

    auto diff = sp_x - sp_direct_x;

    auto error = L2Norm(diff);

    printf("Error: %.2e\n", error);

    // y^t A x
    auto xAy = sp_elem.InnerProduct(x, y);
    auto yAx = sp_elem.InnerProduct(y, x);
    auto symmetric = std::fabs(xAy - yAx) < 1e-8;

    std::cout << std::fabs(xAy - yAx) << std::endl;

    printf("Symmetric: %s\n", symmetric ? "Yes" : "No");

    return EXIT_SUCCESS;
}


//
// This source file is part of appleseed.
// Visit https://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014-2018 Francois Beaune, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// appleseed.foundation headers.
#include "foundation/math/distance.h"
#include "foundation/math/knn.h"
#include "foundation/math/permutation.h"
#include "foundation/math/rng/distribution.h"
#include "foundation/math/rng/mersennetwister.h"
#include "foundation/math/scalar.h"
#include "foundation/math/vector.h"
#include "foundation/platform/defaulttimers.h"
#include "foundation/platform/timers.h"
#include "foundation/utility/iostreamop.h"
#include "foundation/utility/memory.h"
#include "foundation/utility/test.h"
#include "foundation/utility/stopwatch.h"

// Standard headers.
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

using namespace foundation;
using namespace std;

TEST_SUITE(Foundation_Math_Knn_Tree)
{
    TEST_CASE(Empty_GivenDefaultConstructedTree_ReturnsTrue)
    {
        knn::Tree3d tree;

        EXPECT_TRUE(tree.empty());
    }
}

TEST_SUITE(Foundation_Math_Knn_Builder)
{
    TEST_CASE(Build_GivenZeroPoint_BuildsEmptyTree)
    {
        knn::Tree3d tree;

        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(nullptr, 0);

        EXPECT_TRUE(tree.empty());

        EXPECT_EQ(0, tree.m_points.size());
        EXPECT_EQ(0, tree.m_indices.size());

        ASSERT_EQ(1, tree.m_nodes.size());

        EXPECT_TRUE(tree.m_nodes[0].is_leaf());
        EXPECT_EQ(0, tree.m_nodes[0].get_point_count());
        EXPECT_EQ(0, tree.m_nodes[0].get_point_index());
    }

    TEST_CASE(Build_GivenTwoPoints_BuildsCorrectTree)
    {
        const Vector3d Points[] =
        {
            Vector3d(0.0, 0.0, 0.0),
            Vector3d(1.0, 0.0, 0.0)
        };

        knn::Tree3d tree;

        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(Points, 2);

        ASSERT_EQ(2, tree.m_points.size());
        EXPECT_EQ(Points[0], tree.m_points[0]);
        EXPECT_EQ(Points[1], tree.m_points[1]);

        ASSERT_EQ(2, tree.m_indices.size());
        EXPECT_EQ(0, tree.m_indices[0]);
        EXPECT_EQ(1, tree.m_indices[1]);

        ASSERT_EQ(3, tree.m_nodes.size());

        ASSERT_TRUE(tree.m_nodes[0].is_interior());
        EXPECT_EQ(2, tree.m_nodes[0].get_point_count());
        EXPECT_EQ(0, tree.m_nodes[0].get_point_index());

        ASSERT_TRUE(tree.m_nodes[1].is_leaf());
        EXPECT_EQ(1, tree.m_nodes[1].get_point_count());
        EXPECT_EQ(0, tree.m_nodes[1].get_point_index());

        ASSERT_TRUE(tree.m_nodes[2].is_leaf());
        EXPECT_EQ(1, tree.m_nodes[2].get_point_count());
        EXPECT_EQ(1, tree.m_nodes[2].get_point_index());
    }

    TEST_CASE(Build_GivenEightPoints_GeneratesFifteenNodes)
    {
        const size_t PointCount = 8;

        Vector3d points[PointCount];
        for (size_t i = 0; i < PointCount; ++i)
            points[i] = Vector3d(static_cast<double>(PointCount - i - 1), 0.0, 0.0);

        knn::Tree3d tree;

        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(points, PointCount);

        EXPECT_EQ(8 + 4 + 2 + 1, tree.m_nodes.size());
    }

    TEST_CASE(Build_GivenTwoCoincidentPoints_Terminates)
    {
        const size_t PointCount = 2;
        const Vector3d points[PointCount] = { Vector3d(0.0), Vector3d(0.0) };

        knn::Tree3d tree;

        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(points, PointCount);
    }
}

TEST_SUITE(Foundation_Math_Knn_Answer)
{
    TEST_CASE(Size_AfterZeroInsertion_ReturnsZero)
    {
        knn::Answer<double> answer(3);

        EXPECT_EQ(0, answer.size());
    }

    TEST_CASE(Size_AfterOneInsertion_ReturnsOne)
    {
        knn::Answer<double> answer(3);
        answer.array_insert(42, 12.0);

        EXPECT_EQ(1, answer.size());
    }

    TEST_CASE(Empty_AfterZeroInsertion_ReturnsTrue)
    {
        knn::Answer<double> answer(3);

        EXPECT_TRUE(answer.empty());
    }

    TEST_CASE(Empty_AfterOneInsertion_ReturnsFalse)
    {
        knn::Answer<double> answer(3);
        answer.array_insert(42, 12.0);

        EXPECT_FALSE(answer.empty());
    }

    TEST_CASE(Clear_GivenOneItem_EmptiesAnswer)
    {
        knn::Answer<double> answer(3);
        answer.array_insert(42, 12.0);

        answer.clear();

        EXPECT_TRUE(answer.empty());
    }

    TEST_CASE(Sort_GivenFourItemsInSizeFiveAnswer_SortsItems)
    {
        knn::Answer<double> answer(5);
        answer.array_insert(4, 4.0);
        answer.array_insert(3, 3.0);
        answer.array_insert(1, 1.0);
        answer.array_insert(2, 2.0);

        answer.sort();

        EXPECT_EQ(1, answer.get(0).m_index);
        EXPECT_EQ(2, answer.get(1).m_index);
        EXPECT_EQ(3, answer.get(2).m_index);
        EXPECT_EQ(4, answer.get(3).m_index);
    }

    TEST_CASE(MakeHeap_GivenVector_TransformsVectorToHeap)
    {
        knn::Answer<double> answer(5);
        answer.array_insert(5, 5.0);
        answer.array_insert(1, 1.0);
        answer.array_insert(4, 4.0);
        answer.array_insert(3, 3.0);
        answer.array_insert(2, 2.0);

        answer.make_heap();

        for (size_t i = 0; i < answer.size() / 2; ++i)
        {
            const size_t left = 2 * i + 1;
            const size_t right = left + 1;
            EXPECT_LT(answer.get(i).m_square_dist, answer.get(left).m_square_dist);
            EXPECT_LT(answer.get(i).m_square_dist, answer.get(right).m_square_dist);
        }
    }

    TEST_CASE(Insert_GivenCloserItem_KeepsItem)
    {
        knn::Answer<double> answer(3);
        answer.array_insert(1, 1.0);
        answer.array_insert(3, 3.0);
        answer.array_insert(4, 4.0);
        answer.make_heap();

        answer.heap_insert(2, 2.0);

        answer.sort();

        EXPECT_EQ(1, answer.get(0).m_index);
        EXPECT_EQ(2, answer.get(1).m_index);
        EXPECT_EQ(3, answer.get(2).m_index);
    }
}

TEST_SUITE(Foundation_Math_Knn_Query)
{
    TEST_CASE(Run_GivenEightPointsAndQuerySizeFour_ReturnsFourNearestNeighbors)
    {
        const size_t PointCount = 8;
        const size_t AnswerSize = 4;

        Vector3d points[PointCount];
        for (size_t i = 0; i < PointCount; ++i)
            points[i] = Vector3d(static_cast<double>(PointCount - i), 0.0, 0.0);

        knn::Tree3d tree;
        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(points, PointCount);

        knn::Answer<double> answer(AnswerSize);
        knn::Query3d query(tree, answer);
        query.run(Vector3d(4.5, 0.0, 0.0));

        EXPECT_EQ(4, answer.size());
    }

    TEST_CASE(Run_GivenEightPointsAndMaxSearchDistance_ReturnsFourNearestNeighbors)
    {
        const size_t PointCount = 8;
        const size_t AnswerSize = PointCount;
        const double QueryMaxSquareDistance = square(0.5);

        Vector3d points[PointCount];
        for (size_t i = 0; i < PointCount; ++i)
            points[i] = Vector3d(static_cast<double>(PointCount - i), 0.0, 0.0);

        knn::Tree3d tree;
        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(points, PointCount);

        knn::Answer<double> answer(AnswerSize);
        knn::Query3d query(tree, answer);
        query.run(Vector3d(4.75, 0.0, 0.0), QueryMaxSquareDistance);

        EXPECT_EQ(1, answer.size());
    }

    void generate_random_points(
        MersenneTwister&            rng,
        vector<Vector3d>&           points,
        const size_t                count)
    {
        assert(points.empty());

        points.reserve(count);

        for (size_t i = 0; i < count; ++i)
            points.push_back(rand_vector1<Vector3d>(rng));
    }

    TEST_CASE(Run_GivenMaxSearchDistance_ReturnsCorrectResults)
    {
        const size_t PointCount = 1000;
        const size_t QueryCount = 1000;
        const size_t AnswerSize = 100;
        const double QueryMaxSquareDistance = square(0.1);

        MersenneTwister rng;

        vector<Vector3d> points;
        generate_random_points(rng, points, PointCount);

        knn::Tree3d tree;
        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(&points[0], PointCount);

        knn::Answer<double> full_answer(AnswerSize);
        knn::Query3d full_query(tree, full_answer);

        knn::Answer<double> limited_answer(AnswerSize);
        knn::Query3d limited_query(tree, limited_answer);

        for (size_t i = 0; i < QueryCount; ++i)
        {
            const Vector3d q = rand_vector1<Vector3d>(rng);

            full_query.run(q);
            full_answer.sort();

            limited_query.run(q, QueryMaxSquareDistance);
            limited_answer.sort();

            for (size_t j = 0; j < AnswerSize; ++j)
            {
                if (full_answer.get(j).m_square_dist > QueryMaxSquareDistance)
                {
                    EXPECT_EQ(j, limited_answer.size());
                    break;
                }

                EXPECT_EQ(full_answer.get(j).m_index, limited_answer.get(j).m_index);
                EXPECT_EQ(full_answer.get(j).m_square_dist, limited_answer.get(j).m_square_dist);
            }
        }
    }

    struct SortPointByDistancePredicate
    {
        const vector<Vector3d>&     m_points;
        const Vector3d&             m_q;

        SortPointByDistancePredicate(
            const vector<Vector3d>& points,
            const Vector3d&         q)
          : m_points(points)
          , m_q(q)
        {
        }

        bool operator()(const size_t lhs, const size_t rhs) const
        {
            return
                square_distance(m_q, m_points[lhs]) <
                square_distance(m_q, m_points[rhs]);
        }
    };

    void naive_query(
        const vector<Vector3d>&     points,
        const Vector3d&             q,
        vector<size_t>&             indices)
    {
        assert(indices.size() == points.size());

        SortPointByDistancePredicate pred(points, q);
        sort(indices.begin(), indices.end(), pred);
    }

    bool do_results_match_naive_algorithm(
        const vector<Vector3d>&     points,
        const size_t                answer_size,
        const size_t                query_count,
        function<Vector3d()>        make_query_point)
    {
        knn::Tree3d tree;
        knn::Builder3d builder(tree);
        builder.build<DefaultWallclockTimer>(&points[0], points.size());

        knn::Answer<double> answer(answer_size);
        knn::Query3d query(tree, answer);

        vector<size_t> ref_answer(points.size());
        identity_permutation(ref_answer.size(), &ref_answer[0]);

        for (size_t i = 0; i < query_count; ++i)
        {
            const Vector3d q = make_query_point();

            naive_query(points, q, ref_answer);

            query.run(q);
            answer.sort();

            if (answer.size() != answer_size)
                return false;

            for (size_t j = 0; j < answer_size; ++j)
            {
                if (tree.remap(answer.get(j).m_index) != ref_answer[j])
                    return false;
            }
        }

        return true;
    }

    TEST_CASE(Run_UniformPointDistribution_ReturnsIdenticalResultsAsNaiveAlgorithm)
    {
        const size_t PointCount = 1000;
        const size_t QueryCount = 200;
        const size_t AnswerSize = 20;

        MersenneTwister rng;

        vector<Vector3d> points;
        generate_random_points(rng, points, PointCount);

        auto make_query_point = [&rng]() { return rand_vector1<Vector3d>(rng); };
        EXPECT_TRUE(do_results_match_naive_algorithm(points, AnswerSize, QueryCount, make_query_point));
    }

    TEST_CASE(Run_SkewedPointDistribution_ReturnsIdenticalResultsAsNaiveAlgorithm)
    {
        const size_t PointCount = 1000;
        const size_t QueryCount = 200;
        const size_t AnswerSize = 20;

        MersenneTwister rng;

        vector<Vector3d> points;
        generate_random_points(rng, points, PointCount);

        points[0] = Vector3d(0.0);

        for (size_t i = 1; i < points.size(); ++i)
            points[i] = Vector3d(0.55) + 0.45 * points[i];

        auto make_query_point = [&rng]() { return rand_vector1<Vector3d>(rng); };
        EXPECT_TRUE(do_results_match_naive_algorithm(points, AnswerSize, QueryCount, make_query_point));
    }

    TEST_CASE(Run_QueryPointsArePointsFromTheDataSet_ReturnsIdenticalResultsAsNaiveAlgorithm)
    {
        const size_t PointCount = 1000;
        const size_t QueryCount = 200;
        const size_t AnswerSize = 20;

        MersenneTwister rng;

        vector<Vector3d> points;
        generate_random_points(rng, points, PointCount);

        auto make_query_point = [&rng, &points]() { return points[rand_int1(rng, 0, static_cast<int32>(points.size()) - 1)]; };
        EXPECT_TRUE(do_results_match_naive_algorithm(points, AnswerSize, QueryCount, make_query_point));
    }

    struct NeighborData
    {
        size_t                  m_lookup_point;
        vector<size_t>          m_neighbors;
    };

    struct SampleMesh
    {
        vector<Vector3f>        m_points;
        size_t                  m_index_begin;   // inclusive
        size_t                  m_index_end;     // exclusive
        knn::Tree3f             m_tree;
        vector<NeighborData>    m_neighbor_data;
    };

    TEST_CASE(SamTest)
    {
#if 1
        static const char* MiscDataFilePath = "unit tests/inputs/miscData3.txt";
        static const char* SampleDataFilePath = "unit tests/inputs/sampleData3.txt";
        static const char* NeighboursDataFilePath = "unit tests/inputs/neighboursData3.txt";
#elif 0
        static const char* MiscDataFilePath = "unit tests/inputs/miscData2.txt";
        static const char* SampleDataFilePath = "unit tests/inputs/sampleData2.txt";
        static const char* NeighboursDataFilePath = "unit tests/inputs/neighboursData2.txt";
#elif 0
        static const char* MiscDataFilePath = "unit tests/inputs/miscData1.txt";
        static const char* SampleDataFilePath = "unit tests/inputs/sampleData1.txt";
        static const char* NeighboursDataFilePath = "unit tests/inputs/neighboursData1.txt";
#endif

        Stopwatch<DefaultWallclockTimer> sw;

        //
        // Read miscData.txt.
        //

        float radius;

        {
            FILE* f = fopen(MiscDataFilePath, "rt");
            assert(f);

            char token[100];
            fscanf(f, "%s %f\n", token, &radius);
            fclose(f);
        }

        //
        // Read sampleData.txt.
        //

        vector<SampleMesh> sample_meshes;

        {
            FILE* f = fopen(SampleDataFilePath, "rt");
            assert(f);

            SampleMesh current_mesh;
            current_mesh.m_index_begin = 0;
            current_mesh.m_index_end = ~size_t(0);

            size_t current_index = 0;

            while (!feof(f))
            {
                char buf[16 * 1024];
                if (fgets(buf, sizeof(buf), f) == nullptr)
                    break;

                Vector3f p;
                if (sscanf(buf, "%f %f %f\n", &p.x, &p.y, &p.z) == 3)
                {
                    current_mesh.m_points.push_back(p);
                    ++current_index;
                }
                else
                {
                    if (!current_mesh.m_points.empty())
                    {
                        current_mesh.m_points.shrink_to_fit();
                        current_mesh.m_index_end = current_mesh.m_index_begin + current_mesh.m_points.size();
                        sample_meshes.push_back(current_mesh);
                    }

                    current_mesh = SampleMesh();
                    current_mesh.m_index_begin = current_index;
                }
            }

            if (!current_mesh.m_points.empty())
                sample_meshes.push_back(current_mesh);

            fclose(f);
            sample_meshes.shrink_to_fit();
        }

        //
        // Read neighboursData.txt.
        //

        {
            FILE* f = fopen(NeighboursDataFilePath, "rt");
            assert(f);

            while (!feof(f))
            {
                char buf[16 * 1024];
                if (fgets(buf, sizeof(buf), f) == nullptr)
                    break;

                char* ptr = buf;
                int offset = 0;

                NeighborData nd;
                if (sscanf(ptr, "%llu%n", &nd.m_lookup_point, &offset) != 1)
                    break;

                ptr += offset;

                SampleMesh* sample_mesh = &sample_meshes[0];
                while (sample_mesh->m_index_end <= nd.m_lookup_point)
                {
                    ++sample_mesh;
                    assert(sample_mesh < &sample_meshes[0] + sample_meshes.size());
                }

                assert(nd.m_lookup_point >= sample_mesh->m_index_begin);
                assert(nd.m_lookup_point < sample_mesh->m_index_end);

                while (true)
                {
                    size_t index;
                    if (sscanf(ptr, "%llu%n", &index, &offset) != 1)
                        break;

                    assert(index >= sample_mesh->m_index_begin);
                    assert(index < sample_mesh->m_index_end);

                    ptr += offset;
                    nd.m_neighbors.push_back(index - sample_mesh->m_index_begin);
                }

                nd.m_neighbors.shrink_to_fit();
                sort(nd.m_neighbors.begin(), nd.m_neighbors.end());

                sample_mesh->m_neighbor_data.push_back(nd);
            }

            fclose(f);

            for (SampleMesh& sample_mesh : sample_meshes)
                sample_mesh.m_neighbor_data.shrink_to_fit();
        }

        //
        // Build trees.
        //

        sw.start();

        for (SampleMesh& sample_mesh : sample_meshes)
        {
            knn::Builder3f builder(sample_mesh.m_tree);
            builder.build<DefaultWallclockTimer>(&sample_mesh.m_points[0], sample_mesh.m_points.size());
        }

        cerr << "construction in " << sw.measure().get_seconds() * 1000.0 << " ms" << endl;

        //
        // Perform queries.
        //

        const float query_max_square_distance = radius * radius;

        size_t max_answer_size = 0;

        for (const SampleMesh& sample_mesh : sample_meshes)
        {
            for (const NeighborData& nd : sample_mesh.m_neighbor_data)
                max_answer_size = max(max_answer_size, nd.m_neighbors.size());
        }

        // Add one because we will always find the lookup point itself as a neighbor (at distance 0).
        ++max_answer_size;

        max_answer_size = next_pow2(max_answer_size);

        cerr << "max_answer_size = " << max_answer_size << endl;

        size_t query_count = 0, dummy = 0;
        knn::Answer<float> answer(max_answer_size);
        vector<size_t> answer_indices;

        sw.start();

        for (const SampleMesh& sample_mesh : sample_meshes)
        {
            for (const NeighborData& nd : sample_mesh.m_neighbor_data)
            {
                assert(nd.m_lookup_point >= sample_mesh.m_index_begin);
                assert(nd.m_lookup_point < sample_mesh.m_index_end);

                const size_t local_lookup_point_index = nd.m_lookup_point - sample_mesh.m_index_begin;
                const Vector3f& lookup_point = sample_mesh.m_points[local_lookup_point_index];

                knn::Query3f query(sample_mesh.m_tree, answer);
                query.run(lookup_point, query_max_square_distance);
                ++query_count;
                dummy += answer.size();

#ifdef DEBUG
                clear_keep_memory(answer_indices);
                for (size_t i = 0, e = answer.size(); i < e; ++i)
                {
                    const auto& entry = answer.get(i);
                    assert(entry.m_square_dist <= query_max_square_distance);

                    const size_t point_index = sample_mesh.m_tree.remap(entry.m_index);
                    answer_indices.push_back(point_index);

                    const float neighbor_square_dist = square_distance(sample_mesh.m_points[point_index], lookup_point);
                    assert(neighbor_square_dist <= query_max_square_distance);
                }

                const auto self_it = find(answer_indices.begin(), answer_indices.end(), local_lookup_point_index);
                assert(self_it != answer_indices.end());
                answer_indices.erase(self_it);

                sort(answer_indices.begin(), answer_indices.end());

                if (answer_indices.size() != nd.m_neighbors.size())
                {
                    cerr << "query_max_square_distance: " << query_max_square_distance << endl;
                    cerr << "point: " << nd.m_lookup_point << endl;

                    cerr << "answer:" << endl;
                    for (const size_t point_index : answer_indices)
                    {
                        const Vector3f& neighbor_point = sample_mesh.m_points[point_index];
                        const size_t neighbor_index = point_index + sample_mesh.m_index_begin;
                        const float neighbor_square_dist = square_distance(neighbor_point, lookup_point);
                        cerr << "  " << neighbor_index << "\t\t" << neighbor_square_dist << endl;
                    }

                    cerr << "ref:" << endl;
                    for (const size_t point_index : nd.m_neighbors)
                    {
                        const Vector3f& neighbor_point = sample_mesh.m_points[point_index];
                        const size_t neighbor_index = point_index + sample_mesh.m_index_begin;
                        const float neighbor_square_dist = square_distance(neighbor_point, lookup_point);
                        cerr << "  " << neighbor_index << "\t\t" << neighbor_square_dist << endl;
                    }
                }
                else
                {
                    assert(answer_indices == nd.m_neighbors);
                }
#endif
            }
        }

        cerr << query_count << " queries in " << sw.measure().get_seconds() * 1000.0 << " ms (dummy = " << dummy << ")" << endl;
    }
}

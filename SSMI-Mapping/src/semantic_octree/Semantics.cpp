#include <semantic_octree/Semantics.h>
#include <vector>
#include <algorithm>

#define logALPHA static_cast<float>(0)//logALPHA static_cast<float>(log(0.8))
#define log_minus_ALPHA static_cast<float>(0)// log_minus_ALPHA static_cast<float>(log(0.2))

#define logBETA static_cast<float>(0) //static_cast<float>(log(0.1))
#define log_minus_BETA static_cast<float>(0) //log_minus_BETA static_cast<float>(log(0.9))

namespace octomap
{
    // Struct ColorWithLogOdds implementation -------------------------------------
    std::ostream& operator<<(std::ostream& out, ColorWithLogOdds const& c)
    {
        return out << '(' << c.color << ' ' << c.logOdds << ')';
    }

    // Struct SemanticsLogOdds implementation  --------------------------------------
    SemanticsLogOdds SemanticsLogOdds::initSemantics(const ColorOcTreeNode::Color obs, float value,
                                                     float phi, float psi, float maxLogOdds, float minLogOdds)
    {
        SemanticsLogOdds output;
        output.data[0] = ColorWithLogOdds(obs, logBETA + value + psi);
        if(output.data[0].logOdds < minLogOdds)
        {
            output.data[0].logOdds = minLogOdds;
        } else if(output.data[0].logOdds > maxLogOdds)
        {
            output.data[0].logOdds = maxLogOdds;
        }
        
        output.others = log_minus_BETA + value + phi;
        if(output.others < minLogOdds)
        {
            output.others = minLogOdds;
        }
        
        return output;
    }


    SemanticsLogOdds SemanticsLogOdds::semanticFusion(const SemanticsLogOdds l1, const SemanticsLogOdds l2)
    {
        std::vector<ColorOcTreeNode::Color> v1, v2, v3;
        for(int i = 0; i < NUM_SEMANTICS; i++)
        {
            if(l1.data[i].color != ColorOcTreeNode::Color(255,255,255))
                v1.push_back(l1.data[i].color);

            if(l2.data[i].color != ColorOcTreeNode::Color(255,255,255))
                v2.push_back(l2.data[i].color);
        }

        v3 = v1;
        for(int i = 0; i < v2.size(); i++)
        {
            bool exists = false;
            for(int j = 0; j < v1.size(); j++)
            {
                if(v2[i] == v1[j])
                {
                    exists = true;
                    break;
                }
            }
            if(exists == false)
                v3.push_back(v2[i]);
        }

        std::vector<ColorWithLogOdds> avg_data(v3.size());
        float others1 = l1.others - log(v3.size() - v1.size() + 1);
        float others2 = l2.others - log(v3.size() - v2.size() + 1);

        for(int i = 0; i < v3.size(); i++)
        {
            avg_data[i].color = v3[i];
            float logOdds2;
            bool assigned1 = false, assigned2 = false;
            for(int j = 0; j < NUM_SEMANTICS; j++)
            {
                if(avg_data[i].color == l1.data[j].color)
                {
                    avg_data[i].logOdds = l1.data[j].logOdds;
                    assigned1 = true;
                }

                if(avg_data[i].color == l2.data[j].color)
                {
                    logOdds2 = l2.data[j].logOdds;
                    assigned2 = true;
                }
            }

            if(!assigned1)
                avg_data[i].logOdds = others1;

            if(!assigned2)
                logOdds2 = others2;

            avg_data[i].logOdds = (avg_data[i].logOdds + logOdds2) / 2;
        }
        others1 = (others1 + others2) / 2;

        std::sort(avg_data.begin(), avg_data.end()); // Ascending sort
        SemanticsLogOdds output;
        if(avg_data.size() <= NUM_SEMANTICS)
        {
            for(int i = 0; i < avg_data.size(); i++)
                output.data[i] = avg_data[avg_data.size() - 1 - i];

            output.others = others1;
        } else {
            float exp_others = exp(others1);
            for(int i = 0; i < avg_data.size(); i++)
            {
                if(i < NUM_SEMANTICS)
                {
                    output.data[i] = avg_data[avg_data.size() - 1 - i];
                } else {
                    exp_others += exp(avg_data[avg_data.size() - 1 - i].logOdds);
                }
            }

            output.others = log(exp_others);
        }

        return output;
    }

    // nose is semantic and observation is semantic
    SemanticsLogOdds SemanticsLogOdds::fuseObs(const SemanticsLogOdds l, const ColorOcTreeNode::Color obs,
                             float phi, float psi, float maxLogOdds, float minLogOdds)
    {
        std::vector<ColorWithLogOdds> v;
        float others = l.others;
        bool isOthers = true;
        for(int i = 0; i < NUM_SEMANTICS; i++)
        {
            v.push_back(l.data[i]); // jrcv. This line was originaly below
            if(l.data[i].color != ColorOcTreeNode::Color(255,255,255))
            {
                // v.push_back(l.data[i]); 
                // std::cout << "color with logodss:" << i << ":" << l.data[i] << std::endl;
                if(obs == l.data[i].color)
                {
                    isOthers = false;
                    v.back().logOdds += psi;
                    if(v.back().logOdds > maxLogOdds)
                        v.back().logOdds = maxLogOdds;
                } else {
                    v.back().logOdds += phi;
                    if(v.back().logOdds < minLogOdds)
                        v.back().logOdds = minLogOdds;
                }
            }
        }

        if(isOthers)
        {   
            ColorWithLogOdds aux(obs, l.others + logALPHA + psi); // Assigning color and logodss
            v.push_back(aux);
            others += log_minus_ALPHA + phi;

            std::sort(v.begin(), v.end()); //ascending order
            
            if(v.size() > NUM_SEMANTICS)
                {
                    others = log(exp(others) + exp(v[0].logOdds));
                    v.erase(v.begin());
                }
            
        } else {
            others += phi;

            std::sort(v.begin(), v.end()); //ascending order
        }

        if(others < minLogOdds)
        {
            others = minLogOdds;
        } else if(others > maxLogOdds)
        {
            others = maxLogOdds;
        }

        SemanticsLogOdds output;
        output.others = others;
        for(int i = 0; i < v.size(); i++)
            output.data[i] = v[v.size() - 1 - i]; //descending order
        
        return output;
    }
    
    // node is semantic and observation is free
    SemanticsLogOdds SemanticsLogOdds::fuseObsFree(const SemanticsLogOdds l, float phi, float minLogOdds)
    {
        std::vector<ColorWithLogOdds> v;
        SemanticsLogOdds output;
        
        output.others = l.others;
        
        for(int i = 0; i < NUM_SEMANTICS; i++)
        {
            v.push_back(l.data[i]); // jrcv
            
            v.back().logOdds += phi;
            // if(v.back().logOdds < -0.1)
            //     v.back().color = ColorOcTreeNode::Color(255,255,255); //recover free color
            if(v.back().logOdds < minLogOdds)
                v.back().logOdds = minLogOdds;
            // output.data[i].logOdds += phi;
            // if(output.data[i].logOdds < minLogOdds)
            //     output.data[i].logOdds = minLogOdds;
        }
        std::sort(v.begin(), v.end()); //ascending order, jrcv
        
        output.others += phi;
        if(output.others < minLogOdds)
            output.others = minLogOdds;
        for(int i = 0; i < v.size(); i++)
            output.data[i] = v[v.size() - 1 - i]; //descending order
        return output;
    }
    // original implementation of fuseObsFree
    // SemanticsLogOdds SemanticsLogOdds::fuseObsFree(const SemanticsLogOdds l, float phi, float minLogOdds)
    // {
    //     SemanticsLogOdds output = l;
        
    //     for(int i = 0; i < NUM_SEMANTICS; i++)
    //     {
    //         if(output.data[i].color != ColorOcTreeNode::Color(255,255,255))
    //         {
    //             output.data[i].logOdds += phi;
    //             if(output.data[i].logOdds < minLogOdds)
    //                 output.data[i].logOdds = minLogOdds;
    //         }
            
    //     }

    //     output.others += phi;
    //     if(output.others < minLogOdds)
    //         output.others = minLogOdds;
        
    //     return output;
    // }

    std::ostream& operator<<(std::ostream& out, SemanticsLogOdds const& s) {
        out << '(';
        for(int i = 0; i < NUM_SEMANTICS - 1; i++)
            out << s.data[i] << ' ';
        out << s.data[NUM_SEMANTICS - 1];
        out << ')';
        return out;
    }

    /*std::vector<float> softmax(const SemanticsBayesian l)
    {
        std::vector<float> prob_vec (1);
        for(int i = 0; i < NUM_SEMANTICS; i++)
        {
            if(l.data[i].color != ColorOcTreeNode::Color(255,255,255))
                prob_vec.push_back(exp(l.data[i].logOdds));
        }
        prob_vec.push_back(exp(l.others));
        float denom = 0.;
        for (auto& p : prob_vec)
            denom += p;
        for (auto& p : prob_vec)
            p = p / denom;
        return prob_vec;
    }*/

} //namespace octomap

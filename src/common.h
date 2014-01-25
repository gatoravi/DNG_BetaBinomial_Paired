#include<utility>
#include<iostream>
bool compareBaseCounts (std::pair<string, int> i, std::pair<string, int> j) 
{ 
  return (i.second > j.second); 
}

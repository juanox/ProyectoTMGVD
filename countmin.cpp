#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "MurmurHash3.h"
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdint>
#include <limits>
#include <map>
#include<bits/stdc++.h>
#include <utility>
#include <type_traits>
#include <typeinfo>
#include <cxxabi.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main(){
    ifstream inFile;
    inFile.open("ERR626210.fastq");
    string line;
    const int d=4;  // filas de funciones hash.
    const int w=65536;  // # de buckets para cada fila. ///{2^15= 32768 y 2^16=65536}
    const int UMBRAL=10000; // umbral de frecuencia para HH´s. //{CountminCU entre 5000 y 8000, Countmin 16000}
    const int q=40;     ///{q-mer size: 20, 30, 40}
    const int limit=166-q;
    const int DEFAULT=0;
    uint16_t hash_otpt[0];
    vector<vector<int>> C(d, vector<int>(w, DEFAULT));
    unordered_map<string, int> hhitters;
    unordered_map<string, int> count_hhitters;
    unordered_map<string, int>:: iterator itr;
    int cont=0;
    bool flag=false;
    int freq_est;
    // 1era pasada para calcular HH's
    auto start = high_resolution_clock::now();
    while (inFile){ // llenado de matriz de contadores.
        getline(inFile, line);
        if(line[0]=='@'){
            continue;
        }
        else if(line[0]=='+'){
            flag=true;
            continue;
        }
        else if(flag==true){
            flag=false;
            continue;
        }
        else{
            cont++;
            if(line.size()==0){
                break;
            }
            for(int i = 0; i <= limit; i++){
                cout << line << endl;
                string temp = line.substr(i, q);   //  substr(i, q-mer size)
                char qmer[temp.size()+1];
                strcpy(qmer, temp.c_str());
                for(int seed = 0; seed < d; seed++){
                    MurmurHash3_x86_32(qmer, (uint16_t)strlen(qmer), seed, hash_otpt);  //  valor hashed
                    C[seed][hash_otpt[0]]++;    //  se incrementan todos los contadores
                    if(seed == 0){
                        freq_est = C[seed][hash_otpt[0]];
                    }
                    else{
                        if(C[seed][hash_otpt[0]] < freq_est){
                            freq_est = C[seed][hash_otpt[0]];   // se mantiene el minimo en las d filas para el valor hashed
                        }
                    }
                }
                if(freq_est >= UMBRAL){ //  agregar HH´s por sobre UMBRAL = 15000
                    hhitters[qmer] = freq_est;
                }
            }
       }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << '\n' << "Ejecucion: "<< (duration.count())/60 << " minutos, "<< (duration.count())%60 << " segundos.\n";
    cout << "\nLineas procesadas: "<< cont << endl;
    cout << "\n" << "Resultados con CountMin, " << "Parametros: d = "<< d << " w = "<< w << " q = "<< q <<"\n";
    cout << "\n" <<hhitters.size() << " H. Hitters con umbral de " << UMBRAL << "\n\n";

    for (itr = hhitters.begin(); itr != hhitters.end(); itr++){
        cout << itr->first << "  " << itr->second << endl;
    }
    for (itr = hhitters.begin(); itr != hhitters.end(); itr++){ //  seteo inicial de contadores de frecuencias reales en 0's
        count_hhitters[itr->first] = 0;
    }
    // 2da pasada para calculo de frecuencias reales
    inFile.clear();
    inFile.seekg(0,std::ios::beg);
    while (inFile){
        getline(inFile, line);
        if(line[0]=='@'){
            continue;
        }
        else if(line[0]=='+'){
            flag=true;
            continue;
        }
        else if (flag==true){
            flag=false;
            continue;
        }
        else{
            if(line.size()==0){
                break;
            }
            for(int i = 0; i <= limit; i++){
                string temp = line.substr(i, q);   //  substr(i, q-mer size)
                char qmer[temp.size()+1];
                strcpy(qmer, temp.c_str());
                if(count_hhitters.find(qmer) != count_hhitters.end()){
                    count_hhitters.at(qmer)++;
                }
            }
       }
    }
    std::cout << "\nFrecuencias reales de los H. Hitters: \n" << endl;
    for (itr = count_hhitters.begin(); itr != count_hhitters.end(); itr++){
        cout << itr->first << "  " << itr->second << endl;
    }
    inFile.close();
    return 0;
}

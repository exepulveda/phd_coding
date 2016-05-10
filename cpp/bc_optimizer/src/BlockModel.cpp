#include "BlockModel.h"

BlockModel::BlockModel(int _nx,int _ny, int _nz): nx(_nx), ny(_ny), nz(_nz)
{
    this->nxy= nx*ny;
    this->n = nxy*nz;
    this->sizex = 1.0;
    this->sizey = 1.0;
    this->sizez = 1.0;
    this->volume = 1.0;
}

BlockModel::BlockModel()
{
    this->nxy = 0;
    this->n = 0;
    this->sizex = 1.0;
    this->sizey = 1.0;
    this->sizez = 1.0;
    this->volume = 1.0;
}

BlockModel::~BlockModel()
{
    //dtor
}

block_index BlockModel::getBlockId(int i, int j, int k){
    return k*nxy + j*nx + i;
}

double BlockModel::getVolume() {
    return this->volume;
}

void BlockModel::setVolume(float sizex,float sizey,float sizez) {
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->volume = sizex * sizey * sizez;
}

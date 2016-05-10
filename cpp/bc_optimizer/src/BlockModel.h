#ifndef BLOCKMODEL_H
#define BLOCKMODEL_H

typedef unsigned int block_index;

class BlockModel
{
    public:
        BlockModel(int nx,int ny, int nz);
        BlockModel();
        virtual ~BlockModel();
        int nx;
        int ny;
        int nz;
        int nxy;
        int n;
        float sizex;
        float sizey;
        float sizez;

        block_index getBlockId(int i, int j, int k);

        double getVolume();
        void setVolume(float sizex,float sizey,float sizez);

    protected:
    private:
        double volume;

};

#endif // BLOCKMODEL_H

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}


int main()
{
   vtkImageData *image = NewImage(1024, 1350);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   
    //What I added
    int loc, strip;
    for (int i = 0; i < 1350; i++)
    {
        for (int j = 0; j < 1024; j++)
        {
            loc = i*1024*3 + j*3;
            strip = i / 50;
            
            //Blue
            if (strip % 3 == 0)            buffer[loc+2] = 0;
            else if (strip % 3 == 1)       buffer[loc+2] = 128;
            else                           buffer[loc+2] = 255;
            
            //Green
            if ((strip / 3) % 3 == 0)      buffer[loc+1] = 0;
            else if ((strip / 3) % 3 == 1) buffer[loc+1] = 128;
            else                           buffer[loc+1] = 255;
            
            //Red
            if (strip / 9 == 0)            buffer[loc] = 0;
            else if (strip / 9 == 1)       buffer[loc] = 128;
            else                           buffer[loc] = 255;

        }
    }

   WriteImage(image, "proj1A");
}

/***************************************************************************
 *
 *   transformbitmap.cpp -- Rotation and Mirroring of Bitmap
 *
 *   Erik Persson, private, 2000
 */

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include "transformbitmap.h"

HBITMAP TransformBitmap(HBITMAP hBitmap, bool bFlipHorizontal, int nQuartersClockwise)
{
	// Rotate and conditionally flip bitmap
	// Flip is applied before rotation
	
	// Get structure containing bitmap info
	BITMAP info;
	if (GetObject(hBitmap, sizeof(info),&info) != sizeof(info))
		return NULL;
	
	// Get array of data bytes
	int nSrcBytes= info.bmWidthBytes * info.bmHeight;
	BYTE *aSrcBytes= new BYTE[nSrcBytes];
	if (GetBitmapBits(hBitmap, nSrcBytes, aSrcBytes) != nSrcBytes)
	{
		delete[] aSrcBytes;
		return NULL;
	}
	
	// Get source bitmap info
	int bpp = info.bmPlanes * info.bmBitsPixel;		// Bits per pixel
	int nColors= 1 << bpp;
	int nSrcWidth = info.bmWidth;
	int nSrcHeight = info.bmHeight;
	int nSrcRowBytes = ((((nSrcWidth * bpp) + 15) & ~15) / 8);
	
	// Compute dimensions of the resulting bitmap
	int w= nSrcWidth;
	int h= nSrcHeight;
	if (nQuartersClockwise & 1)
	{
		int t= w; w= h; h= t;
	}
	
	// Create an array to hold the result. 
	int nResultRowBytes = ((((w * bpp) + 15) & ~15) / 8);
	// int nResultRowBytes = ((((w * bpp) + 31) & ~31) / 8);
	long nResultBytes = nResultRowBytes * h;
	BYTE *aResultBytes= new BYTE[nResultBytes];
	
	// Set memory to zero
	ZeroMemory( aResultBytes, nResultBytes );
	
	// Now do the actual rotating - a pixel at a time
	// Computing the destination point for each source point
	// will leave a few pixels that do not get covered
	// So we use a reverse transform - e.i. compute the source point
	// for each destination point
	
	// Calculate
	static const src_aX_per_x[4]= {1,0,-1,0};
	static const src_aX_per_y[4]= {0,1,0,-1};
	static const src_aY_per_x[4]= {0,-1,0,1};
	static const src_aY_per_y[4]= {1,-0,-1,0};	
	
	int i= (nQuartersClockwise & 3);
	int src_x_per_x= src_aX_per_x[i];
	int src_x_per_y= src_aX_per_y[i];
	int src_y_per_x= src_aY_per_x[i];
	int src_y_per_y= src_aY_per_y[i];
	
	if (bFlipHorizontal)
	{
		src_x_per_x= -src_x_per_x;
		src_x_per_y= -src_x_per_y;
	}

	// Set starting point
	int src_x= 0, src_y= 0;
	if (src_x_per_x<0 || src_x_per_y<0)
		src_x= nSrcWidth-1; // Start from other end;
	if (src_y_per_x<0 || src_y_per_y<0)
		src_y= nSrcHeight-1; // Start from other end;

	// Compensate for offset added in inner loop
	src_x_per_y -= w*src_x_per_x;
	src_y_per_y -= w*src_y_per_x;
		

	// Each case of bits per pixel is handled separately
	switch(bpp)
	{
	  case 1:		//Monochrome
		  {
			  for (int y= 0; y<h; y++)
			  {
				  for (int x= 0; x<w; x++)
				  {
					  BYTE mask= *(aSrcBytes + nSrcRowBytes*src_y + src_x/8) & (0x80 >> src_x%8);
					  //Adjust mask for destination bitmap
					  mask = mask ? (0x80 >> x%8) : 0;
					  *(aResultBytes + nResultRowBytes*y + x/8) &= ~(0x80 >> x%8);
					  *(aResultBytes + nResultRowBytes*y + x/8) |= mask;
					  src_x += src_x_per_x;
					  src_y += src_y_per_x;
				  }
				  src_x += src_x_per_y;
				  src_y += src_y_per_y;
			  }
		  }
		  break;
	  case 4:
		  {
			  for (int y= 0; y<h; y++)
			  {
				  for (int x= 0; x<w; x++)
				  {
					  BYTE mask= *(aSrcBytes + nSrcRowBytes*src_y + src_x/2) & ((src_x&1) ? 0x0f : 0xf0);
					  //Adjust mask for destination bitmap
					  if( (src_x&1) != (x&1) )
						  mask = (mask&0xf0) ? (mask>>4) : (mask<<4);
					  *(aResultBytes + nResultRowBytes*y + x/2) &= ~((x&1) ? 0x0f : 0xf0);
					  *(aResultBytes + nResultRowBytes*y + x/2) |= mask;
					  src_x += src_x_per_x;
					  src_y += src_y_per_x;
				  }
				  src_x += src_x_per_y;
				  src_y += src_y_per_y;
			  }
		  }
		  break;
	  case 8:
		  {
			  for (int y= 0; y<h; y++)
			  {
				  for (int x= 0; x<w; x++)
				  {
					  aResultBytes[nResultRowBytes*y + x]= 
						  aSrcBytes[nSrcRowBytes*src_y + src_x];
					  
					  src_x += src_x_per_x;
					  src_y += src_y_per_x;
				  }
				  src_x += src_x_per_y;
				  src_y += src_y_per_y;
			  }
		  }
		  break;
	  case 16:
		  {
			  WORD *aSrcWords= (WORD *) aSrcBytes;
			  WORD *aResultWords= (WORD *) aResultBytes;
			  int nSrcRowWords= (nSrcRowBytes>>1);
			  int nResultRowWords= (nResultRowBytes>>1);
			  
			  for (int y= 0; y<h; y++)
			  {
				  for (int x= 0; x<w; x++)
				  {
					  aResultWords[nResultRowBytes*y + x*2]= 
						  aSrcWords[nSrcRowWords*src_y + src_x];
					  
					  src_x += src_x_per_x;
					  src_y += src_y_per_x;
				  }
				  src_x += src_x_per_y;
				  src_y += src_y_per_y;
			  }
		  }
		  break;
	  case 24:
		  {
			  //TRACE("TRANSFORM size (%d,%d): src=(%d,%d)\n",w,h,nSrcWidth,nSrcHeight);
			  for (int y= 0; y<h; y++)
			  {
				  for (int x= 0; x<w; x++)
				  {
					  //TRACE("(%d,%d): src=(%d,%d)\n",x,y,src_x,src_y);
					  DWORD dwPixel= *((LPDWORD)(aSrcBytes + nSrcRowBytes*src_y + src_x*3)) & 0xffffff;
					  *((LPDWORD)(aResultBytes + nResultRowBytes*y + x*3)) |= dwPixel;
					  src_x += src_x_per_x;
					  src_y += src_y_per_x;
				  }
				  src_x += src_x_per_y;
				  src_y += src_y_per_y;
			  }
		  }
		  break;
	  case 32:
		  {
			  DWORD *aSrcDwords= (DWORD *) aSrcBytes;
			  DWORD *aResultDwords= (DWORD *) aResultBytes;
			  int nSrcRowDwords= (nSrcRowBytes>>2);
			  int nResultRowDwords= (nResultRowBytes>>2);
			  
			  for (int y= 0; y<h; y++)
			  {
				  for (int x= 0; x<w; x++)
				  {
					  aResultDwords[nResultRowDwords*y + x]= 
						  aSrcDwords[nSrcRowDwords*src_y + src_x];
					  
					  src_x += src_x_per_x;
					  src_y += src_y_per_x;
				  }
				  src_x += src_x_per_y;
				  src_y += src_y_per_y;
			  }
		  }
		  break;
	}

	// Create a new HBITMAP
	info.bmWidth= w;
	info.bmWidthBytes= nResultRowBytes;
	info.bmHeight= h;
	info.bmBits= aResultBytes;
	HBITMAP hResultBitmap= CreateBitmapIndirect(&info);

	// Free the temporary arrays
	delete[] aSrcBytes;
	delete[] aResultBytes;

	return hResultBitmap;

}

/***************************************************************************
 *
 *   transformbitmap.h -- Rotation and Mirroring of Bitmap
 *
 *   Erik Persson, private, 2000 
 */


#ifndef WINDOWS_TRANSFORMBITMAP_H
#define WINDOWS_TRANSFORMBITMAP_H

HBITMAP TransformBitmap(HBITMAP hBitmap, bool bFlipHorizontal, int nQuartersClockwise);

#endif // WINDOWS_TRANSFORMBITMAP_H
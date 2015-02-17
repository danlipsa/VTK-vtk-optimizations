/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellArray.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCellArray.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkCellArray);

//----------------------------------------------------------------------------
vtkCellArray::vtkCellArray()
{
  this->Points = vtkIdTypeArray::New();
  this->CellIndex = vtkIdTypeArray::New();
  this->CellIndex->InsertValue(0, 0);
  this->InsertPointLocation = 0;
  this->TraversalId = 0;
}

//----------------------------------------------------------------------------
void vtkCellArray::DeepCopy (vtkCellArray *ca)
{
  // Do nothing on a NULL input.
  if (ca == NULL)
    {
    return;
    }

  this->Points->DeepCopy(ca->Points);
  this->CellIndex->DeepCopy(ca->CellIndex);
  this->InsertPointLocation = ca->InsertPointLocation;
  this->TraversalId = ca->TraversalId;
}

//----------------------------------------------------------------------------
vtkCellArray::~vtkCellArray()
{
  this->Points->Delete();
  this->CellIndex->Delete();
}

//----------------------------------------------------------------------------
void vtkCellArray::Initialize()
{
  this->Points->Initialize();
  this->CellIndex->Initialize();
  this->CellIndex->InsertValue(0, 0);
  this->InsertPointLocation = 0;
  this->TraversalId = 0;
}

//----------------------------------------------------------------------------
// Returns the size of the largest cell. The size is the number of points
// defining the cell.
int vtkCellArray::GetMaxCellSize()
{
  int maxSize=0;
  vtkIdType *cellIndex = this->CellIndex->GetPointer(0);
  vtkIdType *maxCellIndex = this->CellIndex->GetPointer(
    this->CellIndex->GetNumberOfTuples() - 1);
  for (; cellIndex < maxCellIndex; ++cellIndex)
    {
    vtkIdType npts;
    if ( (npts= *(cellIndex + 1) - *cellIndex) > maxSize )
      {
      maxSize = npts;
      }
    }
  return maxSize;
}

//----------------------------------------------------------------------------
unsigned long vtkCellArray::GetActualMemorySize()
{
  return this->Points->GetActualMemorySize() + this->CellIndex->GetActualMemorySize();
}

//----------------------------------------------------------------------------
int vtkCellArray::GetNextCell(vtkIdList *pts)
{
  vtkIdType npts, *ppts;
  if (this->GetNextCell(npts, ppts))
    {
    pts->SetNumberOfIds(npts);
    for (vtkIdType i = 0; i < npts; i++)
      {
      pts->SetId(i, ppts[i]);
      }
    return 1;
    }
  return 0;
}

//----------------------------------------------------------------------------
void vtkCellArray::GetCellFromId(vtkIdType cellId, vtkIdList *pts)
{
  vtkIdType npts;
  vtkIdType *ppts;
  this->GetCellFromId(cellId, npts, ppts);
  pts->SetNumberOfIds(npts);
  for (vtkIdType i = 0; i < npts; i++)
    {
    pts->SetId(i, ppts[i]);
    }
}

//----------------------------------------------------------------------------
void vtkCellArray::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Number Of Cells: " << this->GetNumberOfCells() << endl;
  os << indent << "Insert Point Location: " << this->InsertPointLocation << endl;
  os << indent << "Traversal CellId: " << this->TraversalId << endl;
}

//----------------------------------------------------------------------------
int vtkCellArray::IsHomogeneous(int* numberOfPoints) const
{
  int npts = 0;
  vtkIdType *cellIndex = this->CellIndex->GetPointer(0);
  vtkIdType *maxCellIndex = this->CellIndex->GetPointer(
    this->CellIndex->GetNumberOfTuples() - 1);
  if (cellIndex >= maxCellIndex)
    {
    // there are 0 cells, so numberOfPoints is undefined
    return 1;
    }
  npts = *(cellIndex + 1) - *cellIndex;
  ++cellIndex;
  for (; cellIndex < maxCellIndex; ++cellIndex)
    {
    if ((*(cellIndex + 1) - *cellIndex) != npts)
      {
      return 0;
      }
    }
  *numberOfPoints = npts;
  return 1;
}

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXMLUnstructuredDataWriter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkXMLUnstructuredDataWriter.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataCompressor.h"
#include "vtkDataSetAttributes.h"
#include "vtkErrorCode.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"

vtkCxxRevisionMacro(vtkXMLUnstructuredDataWriter, "1.12");

//----------------------------------------------------------------------------
vtkXMLUnstructuredDataWriter::vtkXMLUnstructuredDataWriter()
{
  this->NumberOfPieces = 1;
  this->WritePiece = -1;
  this->GhostLevel = 0;
  this->CellPoints = vtkIdTypeArray::New();
  this->CellOffsets = vtkIdTypeArray::New();
  this->CellPoints->SetName("connectivity");
  this->CellOffsets->SetName("offsets");

  this->CurrentPiece = 0;
}

//----------------------------------------------------------------------------
vtkXMLUnstructuredDataWriter::~vtkXMLUnstructuredDataWriter()
{
  this->CellPoints->Delete();
  this->CellOffsets->Delete();
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "NumberOfPieces: " << this->NumberOfPieces << "\n";
  os << indent << "WritePiece: " << this->WritePiece << "\n";
  os << indent << "GhostLevel: " << this->GhostLevel << "\n";
}

//----------------------------------------------------------------------------
vtkPointSet* vtkXMLUnstructuredDataWriter::GetInputAsPointSet()
{
  return static_cast<vtkPointSet*>(this->Superclass::GetInput());
}

//----------------------------------------------------------------------------
int vtkXMLUnstructuredDataWriter::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{

  if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    if((this->WritePiece < 0) || (this->WritePiece >= this->NumberOfPieces))
      {
      this->SetInputUpdateExtent(
        this->CurrentPiece, this->NumberOfPieces, this->GhostLevel);
      }
    else
      {
      this->SetInputUpdateExtent(
        this->WritePiece, this->NumberOfPieces, this->GhostLevel);
      }
    return 1;
    }
  
  // generate the data
  else if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    this->SetErrorCode(vtkErrorCode::NoError);

    if(!this->Stream && !this->FileName)
      {
      this->SetErrorCode(vtkErrorCode::NoFileNameError);
      vtkErrorMacro("The FileName or Stream must be set first.");
      return 0;
      }

    // We don't want to write more pieces than the pipeline can produce,
    // but we need to preserve the user's requested number of pieces in
    // case the input changes later.  If MaximumNumberOfPieces is lower
    // than 1, any number of pieces can be produced by the pipeline.
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    int numPieces = this->NumberOfPieces;
    int maxPieces = 
      inInfo->Get(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES());
    if((maxPieces > 0) && (this->NumberOfPieces > maxPieces))
      {
      this->NumberOfPieces = maxPieces;
      }

    if (this->WritePiece >= 0)
      {
      this->CurrentPiece = this->WritePiece;
      }
    else
      {
      float progressRange[2] = {0,0};
      this->GetProgressRange(progressRange);
      this->SetProgressRange(progressRange, this->CurrentPiece, this->NumberOfPieces);
      }

    int result = 1;
    if (this->CurrentPiece == 0 || this->WritePiece >= 0)
      {
      // We are just starting to write.  Do not call
      // UpdateProgressDiscrete because we want a 0 progress callback the
      // first time.
      this->UpdateProgress(0);

      // Initialize progress range to entire 0..1 range.
      float wholeProgressRange[2] = {0,1};
      this->SetProgressRange(wholeProgressRange, 0, 1);
  
      if (!this->OpenFile())
        {
        this->NumberOfPieces = numPieces;
        return 0;
        }
      // Write the file.
      if (!this->StartFile())
        {
        this->NumberOfPieces = numPieces;
        return 0;
        }

      if (!this->WriteHeader())
        {
        this->NumberOfPieces = numPieces;
        return 0;
        }

      if( this->DataMode == vtkXMLWriter::Appended && this->FieldDataOffsets)
        {
        // Write the field data arrays.
        this->WriteFieldDataAppendedData(this->GetInput()->GetFieldData(),
          this->FieldDataOffsets);
        if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
          {
          this->DeletePositionArrays();
          return 0;
          }
        }
      }

    result = this->WriteAPiece();

    if((this->WritePiece < 0) || (this->WritePiece >= this->NumberOfPieces))
      {
      // Tell the pipeline to start looping.
      if (this->CurrentPiece == 0)
        {
        request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
        }
      this->CurrentPiece++;
      }

    if (this->CurrentPiece == this->NumberOfPieces || this->WritePiece >= 0)
      {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentPiece = 0;
  
      if (!this->WriteFooter())
        {
        this->NumberOfPieces = numPieces;
        return 0;
        }

      if (!this->EndFile())
        {
        this->NumberOfPieces = numPieces;
        return 0;
        }

      this->CloseFile();
      }
    this->NumberOfPieces = numPieces;

    // We have finished writing.
    this->UpdateProgressDiscrete(1);
    return result;
    }
  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::AllocatePositionArrays()
{
  this->NumberOfPointsPositions = new unsigned long[this->NumberOfPieces];
  this->PointsPositions = new unsigned long[this->NumberOfPieces];
  this->PointDataPositions = new unsigned long*[this->NumberOfPieces];
  this->CellDataPositions = new unsigned long*[this->NumberOfPieces];
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::DeletePositionArrays()
{
  delete [] this->NumberOfPointsPositions;
  delete [] this->PointsPositions;
  delete [] this->PointDataPositions;
  delete [] this->CellDataPositions;
  if (this->FieldDataOffsets)
    {
    delete [] this->FieldDataOffsets;
    this->FieldDataOffsets = NULL;
    }
}


//----------------------------------------------------------------------------
int vtkXMLUnstructuredDataWriter::WriteHeader()
{
  vtkIndent indent = vtkIndent().GetNextIndent();

  ostream& os = *(this->Stream);
  
  // Open the primary element.
  os << indent << "<" << this->GetDataSetName() << ">\n";

  this->WriteFieldData(indent.GetNextIndent());

  if(this->DataMode == vtkXMLWriter::Appended)
    {
    vtkIndent nextIndent = indent.GetNextIndent();

    this->AllocatePositionArrays();

    if((this->WritePiece < 0) || (this->WritePiece >= this->NumberOfPieces))
      {
      // Loop over each piece and write its structure.
      int i;
      for(i=0; i < this->NumberOfPieces; ++i)
        {
        // Open the piece's element.
        os << nextIndent << "<Piece";
        this->WriteAppendedPieceAttributes(i);
        if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
          {
          this->DeletePositionArrays();
          return 0;
          }
        os << ">\n";
        
        this->WriteAppendedPiece(i, nextIndent.GetNextIndent());
        if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
          {
          this->DeletePositionArrays();
          return 0;
          }
        
        // Close the piece's element.
        os << nextIndent << "</Piece>\n";
        }
      }
    else
      {
      // Write just the requested piece.
      // Open the piece's element.
      os << nextIndent << "<Piece";
      this->WriteAppendedPieceAttributes(this->WritePiece);
      if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
        {
        this->DeletePositionArrays();
        return 0;
        }
      os << ">\n";
      
      this->WriteAppendedPiece(this->WritePiece, nextIndent.GetNextIndent());
      if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
        {
        this->DeletePositionArrays();
        return 0;
        }
      
      // Close the piece's element.
      os << nextIndent << "</Piece>\n";
      }

    // Close the primary element.
    os << indent << "</" << this->GetDataSetName() << ">\n";
    os.flush();
    if (os.fail())
      {
      this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
      this->DeletePositionArrays();
      return 0;
      }
    
    this->StartAppendedData();
    if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
      {
      this->DeletePositionArrays();
      return 0;
      }
    
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkXMLUnstructuredDataWriter::WriteAPiece()
{
  vtkIndent indent = vtkIndent().GetNextIndent();

  int result=1;

  if(this->DataMode == vtkXMLWriter::Appended)
    {
    this->WriteAppendedPieceData(this->CurrentPiece);
    }
  else
    {
    result = this->WriteInlineMode(indent);
    }

  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    this->DeletePositionArrays();
    result = 0;
    }
  return result;
}

//----------------------------------------------------------------------------
int vtkXMLUnstructuredDataWriter::WriteFooter()
{
  vtkIndent indent = vtkIndent().GetNextIndent();

  ostream& os = *(this->Stream);
  
  if(this->DataMode == vtkXMLWriter::Appended)
    {
    this->DeletePositionArrays();
    this->EndAppendedData();
    }
  else
    {
    // Close the primary element.
    os << indent << "</" << this->GetDataSetName() << ">\n";
    os.flush();
    if (os.fail())
      {
      return 0;
      }
    }
  
  return 1;
}

//----------------------------------------------------------------------------
int vtkXMLUnstructuredDataWriter::WriteInlineMode(vtkIndent indent)
{
  ostream& os = *(this->Stream);
  vtkIndent nextIndent = indent.GetNextIndent();
  
  // Open the piece's element.
  os << nextIndent << "<Piece";
  this->WriteInlinePieceAttributes();
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return 0;
    }
  os << ">\n";
  
  this->WriteInlinePiece(nextIndent.GetNextIndent());
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return 0;
    }

  // Close the piece's element.
  os << nextIndent << "</Piece>\n";

  return 1;
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::WriteInlinePieceAttributes()
{
  vtkPointSet* input = this->GetInputAsPointSet();
  this->WriteScalarAttribute("NumberOfPoints",
                             input->GetNumberOfPoints());
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::WriteInlinePiece(vtkIndent indent)
{
  vtkPointSet* input = this->GetInputAsPointSet();
  
  // Split progress among point data, cell data, and point arrays.
  float progressRange[2] = {0,0};
  this->GetProgressRange(progressRange);
  float fractions[4];
  this->CalculateDataFractions(fractions);
  
  // Set the range of progress for the point data arrays.
  this->SetProgressRange(progressRange, 0, fractions);  
  
  // Write the point data arrays.
  this->WritePointDataInline(input->GetPointData(), indent);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  // Set the range of progress for the cell data arrays.
  this->SetProgressRange(progressRange, 1, fractions);
  
  // Write the cell data arrays.
  this->WriteCellDataInline(input->GetCellData(), indent);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  // Set the range of progress for the point specification array.
  this->SetProgressRange(progressRange, 2, fractions);
  
  // Write the point specification array.
  this->WritePointsInline(input->GetPoints(), indent);
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::WriteAppendedPieceAttributes(int index)
{
  this->NumberOfPointsPositions[index] =
    this->ReserveAttributeSpace("NumberOfPoints");
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::WriteAppendedPiece(int index,
                                                      vtkIndent indent)
{
  vtkPointSet* input = this->GetInputAsPointSet();
  
  this->PointDataPositions[index] =
    this->WritePointDataAppended(input->GetPointData(), indent);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  this->CellDataPositions[index] =
    this->WriteCellDataAppended(input->GetCellData(), indent);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  this->PointsPositions[index] =
    this->WritePointsAppended(input->GetPoints(), indent);
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::WriteAppendedPieceData(int index)
{
  ostream& os = *(this->Stream);
  vtkPointSet* input = this->GetInputAsPointSet();
  
  unsigned long returnPosition = os.tellp();
  os.seekp(this->NumberOfPointsPositions[index]);
  vtkPoints* points = input->GetPoints();
  this->WriteScalarAttribute("NumberOfPoints",
                             (points?points->GetNumberOfPoints():0));
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  os.seekp(returnPosition);
  
  // Split progress among point data, cell data, and point arrays.
  float progressRange[2] = {0,0};
  this->GetProgressRange(progressRange);
  float fractions[4];
  this->CalculateDataFractions(fractions);
  
  // Set the range of progress for the point data arrays.
  this->SetProgressRange(progressRange, 0, fractions);
  
  // Write the point data arrays.
  this->WritePointDataAppendedData(input->GetPointData(),
                                   this->PointDataPositions[index]);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  // Set the range of progress for the cell data arrays.
  this->SetProgressRange(progressRange, 1, fractions);
  
  // Write the cell data arrays.  
  this->WriteCellDataAppendedData(input->GetCellData(),
                                  this->CellDataPositions[index]);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  // Set the range of progress for the point specification array.
  this->SetProgressRange(progressRange, 2, fractions);
  
  // Write the point specification array.
  this->WritePointsAppendedData(input->GetPoints(),
                                this->PointsPositions[index]);
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::WriteCellsInline(const char* name,
                                                    vtkCellArray* cells,
                                                    vtkDataArray* types,
                                                    vtkIndent indent)
{
  this->ConvertCells(cells);
  
  ostream& os = *(this->Stream);
  os << indent << "<" << name << ">\n";
  
  // Split progress by cell connectivity, offset, and type arrays.
  float progressRange[2] = {0,0};
  this->GetProgressRange(progressRange);
  float fractions[4];
  this->CalculateCellFractions(fractions, types?types->GetNumberOfTuples():0);
  
  // Set the range of progress for the connectivity array.
  this->SetProgressRange(progressRange, 0, fractions);
  
  // Write the connectivity array.
  this->WriteDataArrayInline(this->CellPoints, indent.GetNextIndent());
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  // Set the range of progress for the offsets array.
  this->SetProgressRange(progressRange, 1, fractions);
  
  // Write the offsets array.
  this->WriteDataArrayInline(this->CellOffsets, indent.GetNextIndent());
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    return;
    }
  
  if(types)
    {
    // Set the range of progress for the types array.
    this->SetProgressRange(progressRange, 2, fractions);
    
    // Write the types array.
    this->WriteDataArrayInline(types, indent.GetNextIndent(), "types");
    if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
      {
      return;
      }
    }
  os << indent << "</" << name << ">\n";
  os.flush();
  if (os.fail())
    {
    this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
    }
}

//----------------------------------------------------------------------------
unsigned long*
vtkXMLUnstructuredDataWriter::WriteCellsAppended(const char* name,
                                                 vtkDataArray* types,
                                                 vtkIndent indent)
{
  unsigned long* positions = new unsigned long[3];
  ostream& os = *(this->Stream);
  os << indent << "<" << name << ">\n";
  positions[0] = this->WriteDataArrayAppended(this->CellPoints,
                                              indent.GetNextIndent());
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    delete [] positions;
    return NULL;
    }
  
  positions[1] = this->WriteDataArrayAppended(this->CellOffsets,
                                              indent.GetNextIndent());
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    delete [] positions;
    return NULL;
    }
  
  if(types)
    {
    positions[2] = this->WriteDataArrayAppended(types, indent.GetNextIndent(),
                                                "types");
    if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
      {
      delete [] positions;
      return NULL;
      }
    }
  os << indent << "</" << name << ">\n";
  os.flush();
  if (os.fail())
    {
    this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
    delete [] positions;
    return NULL;
    }
  
  return positions;
}

//----------------------------------------------------------------------------
void
vtkXMLUnstructuredDataWriter::WriteCellsAppendedData(vtkCellArray* cells,
                                                     vtkDataArray* types,
                                                     unsigned long* positions)
{
  this->ConvertCells(cells);
  
  // Split progress by cell connectivity, offset, and type arrays.
  float progressRange[2] = {0,0};
  this->GetProgressRange(progressRange);
  float fractions[4];
  this->CalculateCellFractions(fractions, types?types->GetNumberOfTuples():0);
  
  // Set the range of progress for the connectivity array.
  this->SetProgressRange(progressRange, 0, fractions);
  
  // Write the connectivity array.
  this->WriteDataArrayAppendedData(this->CellPoints, positions[0]);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    delete [] positions;
    return;
    }
  
  // Set the range of progress for the offsets array.
  this->SetProgressRange(progressRange, 1, fractions);
  
  // Write the offsets array.
  this->WriteDataArrayAppendedData(this->CellOffsets, positions[1]);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
    {
    delete [] positions;
    return;
    }
  
  if(types)
    {
    // Set the range of progress for the types array.
    this->SetProgressRange(progressRange, 2, fractions);
    
    // Write the types array.
    this->WriteDataArrayAppendedData(types, positions[2]);
    if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
      {
      delete [] positions;
      return;
      }
    }
  delete [] positions;
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::ConvertCells(vtkCellArray* cells)
{
  vtkIdTypeArray* connectivity = cells->GetData();
  vtkIdType numberOfCells = cells->GetNumberOfCells();
  vtkIdType numberOfTuples = connectivity->GetNumberOfTuples();
  
  this->CellPoints->SetNumberOfTuples(numberOfTuples - numberOfCells);
  this->CellOffsets->SetNumberOfTuples(numberOfCells);
  
  vtkIdType* inCell = connectivity->GetPointer(0);
  vtkIdType* outCellPointsBase = this->CellPoints->GetPointer(0);
  vtkIdType* outCellPoints = outCellPointsBase;
  vtkIdType* outCellOffset = this->CellOffsets->GetPointer(0);
  
  vtkIdType i;
  for(i=0;i < numberOfCells; ++i)
    {
    vtkIdType numberOfPoints = *inCell++;
    memcpy(outCellPoints, inCell, sizeof(vtkIdType)*numberOfPoints);
    outCellPoints += numberOfPoints;
    inCell += numberOfPoints;
    *outCellOffset++ = outCellPoints - outCellPointsBase;
    }
}

//----------------------------------------------------------------------------
vtkIdType vtkXMLUnstructuredDataWriter::GetNumberOfInputPoints()
{
  vtkPointSet* input = this->GetInputAsPointSet();
  vtkPoints* points = input->GetPoints();
  return points?points->GetNumberOfPoints():0;
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::CalculateDataFractions(float* fractions)
{
  // Calculate the fraction of point/cell data and point
  // specifications contributed by each component.
  vtkPointSet* input = this->GetInputAsPointSet();
  int pdArrays = input->GetPointData()->GetNumberOfArrays();
  int cdArrays = input->GetCellData()->GetNumberOfArrays();
  vtkIdType pdSize = pdArrays*this->GetNumberOfInputPoints();
  vtkIdType cdSize = cdArrays*this->GetNumberOfInputCells();
  int total = (pdSize+cdSize+this->GetNumberOfInputPoints());
  if(total == 0)
    {
    total = 1;
    }
  fractions[0] = 0;
  fractions[1] = float(pdSize)/total;
  fractions[2] = float(pdSize+cdSize)/total;
  fractions[3] = 1;
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::CalculateCellFractions(float* fractions,
                                                          vtkIdType typesSize)
{
  // Calculate the fraction of cell specification data contributed by
  // each of the connectivity, offset, and type arrays.
  vtkIdType connectSize = this->CellPoints->GetNumberOfTuples();
  vtkIdType offsetSize = this->CellOffsets->GetNumberOfTuples();
  vtkIdType total = connectSize+offsetSize+typesSize;
  if(total == 0)
    {
    total = 1;
    }
  fractions[0] = 0;
  fractions[1] = float(connectSize)/total;
  fractions[2] = float(connectSize+offsetSize)/total;
  fractions[3] = 1;
}

//----------------------------------------------------------------------------
void vtkXMLUnstructuredDataWriter::SetInputUpdateExtent(
  int piece, int numPieces, int ghostLevel)
{
  vtkInformation* inInfo = 
    this->GetExecutive()->GetInputInformation(0, 0);
  inInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  inInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), ghostLevel);
}


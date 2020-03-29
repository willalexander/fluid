#include <Windows.h>
#include <mfapi.h>
#include <mfidl.h>
#include <Mfreadwrite.h>
#include <mferror.h>

#pragma comment(lib, "mfreadwrite")
#pragma comment(lib, "mfplat")
#pragma comment(lib, "mfuuid")

#include <iostream>


template <class T> void SafeRelease(T **ppT)
{
	if (*ppT)
	{
		(*ppT)->Release();
		*ppT = NULL;
	}
}

// Format constants
UINT32 VIDEO_WIDTH = 400;
UINT32 VIDEO_HEIGHT = 400;
const UINT32 VIDEO_FPS = 25;
const UINT64 VIDEO_FRAME_DURATION = 10 * 1000 * 1000 / VIDEO_FPS;
const UINT32 VIDEO_BIT_RATE = 3200000;
const GUID   VIDEO_ENCODING_FORMAT = MFVideoFormat_WMV3;
const GUID   VIDEO_INPUT_FORMAT = MFVideoFormat_RGB32;
UINT32 VIDEO_PELS = VIDEO_WIDTH * VIDEO_HEIGHT;
UINT32 VIDEO_FRAME_COUNT = 4 * VIDEO_FPS;

// Buffer to hold the video frame data.
DWORD *videoFrameBuffer;

IMFSinkWriter *pSinkWriter = NULL;
DWORD stream;
LONGLONG rtStart = 0;


void videoOut_initialize(char **, int, int, float);
void videoOut_addFrame();
void videoOut_finalize();


HRESULT InitializeSinkWriter(IMFSinkWriter **ppWriter, DWORD *pStreamIndex)
{
	*ppWriter = NULL;
	*pStreamIndex = NULL;

	IMFSinkWriter   *pSinkWriter = NULL;
	IMFMediaType    *pMediaTypeOut = NULL;
	IMFMediaType    *pMediaTypeIn = NULL;
	DWORD           streamIndex;

	HRESULT hr = MFCreateSinkWriterFromURL(L"output.wmv", NULL, NULL, &pSinkWriter);

	// Set the output media type.
	if (SUCCEEDED(hr)) hr = MFCreateMediaType(&pMediaTypeOut);
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetGUID(MF_MT_SUBTYPE, VIDEO_ENCODING_FORMAT);
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetUINT32(MF_MT_AVG_BITRATE, VIDEO_BIT_RATE);
	if (SUCCEEDED(hr)) hr = pMediaTypeOut->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);
	if (SUCCEEDED(hr)) hr = MFSetAttributeSize(pMediaTypeOut, MF_MT_FRAME_SIZE, VIDEO_WIDTH, VIDEO_HEIGHT);
	if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeOut, MF_MT_FRAME_RATE, VIDEO_FPS, 1);
	if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeOut, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
	if (SUCCEEDED(hr)) hr = pSinkWriter->AddStream(pMediaTypeOut, &streamIndex);

	// Set the input media type.
	if (SUCCEEDED(hr)) hr = MFCreateMediaType(&pMediaTypeIn);
	if (SUCCEEDED(hr)) hr = pMediaTypeIn->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
	if (SUCCEEDED(hr)) hr = pMediaTypeIn->SetGUID(MF_MT_SUBTYPE, VIDEO_INPUT_FORMAT);
	if (SUCCEEDED(hr)) hr = pMediaTypeIn->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);
	if (SUCCEEDED(hr)) hr = MFSetAttributeSize(pMediaTypeIn, MF_MT_FRAME_SIZE, VIDEO_WIDTH, VIDEO_HEIGHT);
	if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeIn, MF_MT_FRAME_RATE, VIDEO_FPS, 1);
	if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(pMediaTypeIn, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
	if (SUCCEEDED(hr)) hr = pSinkWriter->SetInputMediaType(streamIndex, pMediaTypeIn, NULL);

	// Tell the sink writer to start accepting data.
	if (SUCCEEDED(hr)) hr = pSinkWriter->BeginWriting();

	// Return the pointer to the caller.
	if (SUCCEEDED(hr))
	{
		*ppWriter = pSinkWriter;
		(*ppWriter)->AddRef();
		*pStreamIndex = streamIndex;
	}

	SafeRelease(&pSinkWriter);
	SafeRelease(&pMediaTypeOut);
	SafeRelease(&pMediaTypeIn);
	return hr;
}

HRESULT WriteFrame(
	IMFSinkWriter *pWriter,
	DWORD streamIndex,
	const LONGLONG& rtStart        // Time stamp.
)
{
	IMFSample *pSample = NULL;
	IMFMediaBuffer *pBuffer = NULL;

	const LONG cbWidth = 4 * VIDEO_WIDTH;
	const DWORD cbBuffer = cbWidth * VIDEO_HEIGHT;

	BYTE *pData = NULL;

	// Create a new memory buffer.
	HRESULT hr = MFCreateMemoryBuffer(cbBuffer, &pBuffer);

	// Lock the buffer and copy the video frame to the buffer.
	if (SUCCEEDED(hr)) hr = pBuffer->Lock(&pData, NULL, NULL);
	if (SUCCEEDED(hr))
	{
		hr = MFCopyImage(
			pData,                      // Destination buffer.
			cbWidth,                    // Destination stride.
			(BYTE*)videoFrameBuffer,    // First row in source image.
			cbWidth,                    // Source stride.
			cbWidth,                    // Image width in bytes.
			VIDEO_HEIGHT                // Image height in pixels.
		);
	}
	if (pBuffer)
	{
		pBuffer->Unlock();
	}

	// Set the data length of the buffer.
	if (SUCCEEDED(hr)) hr = pBuffer->SetCurrentLength(cbBuffer);

	// Create a media sample and add the buffer to the sample.
	if (SUCCEEDED(hr)) hr = MFCreateSample(&pSample);
	if (SUCCEEDED(hr)) hr = pSample->AddBuffer(pBuffer);

	// Set the time stamp and the duration.
	if (SUCCEEDED(hr)) hr = pSample->SetSampleTime(rtStart);
	if (SUCCEEDED(hr)) hr = pSample->SetSampleDuration(VIDEO_FRAME_DURATION);

	// Send the sample to the Sink Writer.
	if (SUCCEEDED(hr)) hr = pWriter->WriteSample(streamIndex, pSample);

	SafeRelease(&pSample);
	SafeRelease(&pBuffer);
	return hr;
}

void videoOut_initialize(char **hexBufferOut, int vw, int vh, float d)
{
	VIDEO_WIDTH = (UINT32)(vw);
	VIDEO_HEIGHT = (UINT32)(vh);
	VIDEO_FRAME_COUNT = (UINT32)(d * 25);
	VIDEO_PELS = VIDEO_WIDTH * VIDEO_HEIGHT;

	videoFrameBuffer = (DWORD *)(malloc(VIDEO_PELS * sizeof(DWORD)));

	HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
	if (SUCCEEDED(hr))
	{
		hr = MFStartup(MF_VERSION);
		if (SUCCEEDED(hr))
		{
			hr = InitializeSinkWriter(&pSinkWriter, &stream);
			rtStart = 0;
		}
	}

	*hexBufferOut = (char *)(videoFrameBuffer);
}

void videoOut_addFrame()
{
	HRESULT hr = WriteFrame(pSinkWriter, stream, rtStart);

	rtStart += VIDEO_FRAME_DURATION;
}

void videoOut_finalize()
{
	HRESULT hr = pSinkWriter->Finalize();
	SafeRelease(&pSinkWriter);
	MFShutdown();
	CoUninitialize();
	free(videoFrameBuffer);
}
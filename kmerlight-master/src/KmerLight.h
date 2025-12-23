/*
 * KmerLight.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_KMERLIGHT_H_
#define SRC_KMERLIGHT_H_
#include <pthread.h>
#include <time.h>
#include "types.h"
#include "Kmer.h"
#include "constants.h"
#include "FileRead.h"
#include "KmerHash.h"
#include "CountSketch.h"

#define nbufs 2 

struct KMLThrdData{
	byte_t * buf;
	int hashbufidx;
	int start;
	int end;
	unsigned int size;
	int id;
	CountSketch * cs = NULL;
};


class KmerLight{
	static byte_t * seqbufs[nbufs]; //存储从文件中读取的数据

	static Rands rands[Constants::copies];//随机生成器，用于随机化


public:


	CountSketch CS;//用于估计k-mer的丰度

	static void (*partAdjust) (byte_t * & seqbuf, int & start, int & end, int & idx );
	static void (*bufAdjust) (byte_t * & seqbuf,  int & start,  int end);
     static int (*streamSkip) (byte_t * & seqbuf, int & start, int & end, int & idx); //个函数用于在处理数据时调整缓冲区和跳过序列的注释行


	static FileRead filep; //读取文件

	static pthread_t thrds[3];
	static KMLThrdData thrd_data[3]; //定义了一组线程和线程数据

	static unsigned long memsize;  
	static int readBeg;  //用于跟踪内存使用情况和读取位置




	KmerLight(int k,  char ** infiles, int nfiles){
		KmerHash::init(k);
		for(int i=0; i< nbufs; ++i) {//遍历缓冲区的数量，为每个缓冲区分配内存
			seqbufs[i] = (byte_t *) new byte_t[Constants::seq_buf_size];//char的数组定义
			memsize += Constants::seq_buf_size;//用于跟踪使用的总内存量
		}
		KmerLight::memsize += (1+Constants::buckets)*4*Constants::R*Constants::copies;//R是桶大小，bucket是桶的数量

		KmerLight::memsize += (16*Constants::seq_buf_size * KmerHash::km_word_len);
		filep.reset();
		filep.openFiles(infiles, nfiles);

		switch(filep.getType()){
			case FASTA:
				partAdjust = partAdjustFasta;
				bufAdjust = bufAdjustFasta;
				streamSkip = streamSkipFasta;
				break;
			case FASTQ:
				partAdjust = partAdjustFastQ;
				streamSkip = streamSkipFastaQ;
				bufAdjust = bufAdjustFastQ;
				break;
		}
	}
	~KmerLight(){
		for(int i=0; i< nbufs; ++i) {
			free (seqbufs[i]);
		}
	}

	static size_t loadStream(byte_t* seqbuf){//从文件加载数据到缓冲区
		return filep.fileRead(seqbuf);
	}

	static void *loadStreamWrap(void *inp){//从文件中加载数据
		KMLThrdData * data = (KMLThrdData *)inp;//将void参数转换成KMLThrdData类型
		unsigned int size;
		size = (unsigned int) loadStream(data->buf);//调用loadStream函数，
		data->size = size;
		data->start = 0;
		pthread_exit(NULL);
	}

	static void *storeValsWrap(void *inp){//线程函数
		KMLThrdData * data = (KMLThrdData *)inp;
		(*data->cs).storeValsInSketchMT(data->hashbufidx);//存储哈希之后的值了
		pthread_exit(NULL);
	}

	static void *streamHashPartWrap(void *inp){//用于处理序列数据的一部分
			KMLThrdData * data = (KMLThrdData *)inp;
			 int start = data -> start;
			int size = data -> end;
			bufAdjust(data->buf, start, size);//根据需要处理的部分调整缓冲区
			data-> start = start;
			data-> end = size;
			streamHashPart(data->buf,data->start, data->end, data->hashbufidx);//存储k-mer
			pthread_exit(NULL);
	}
 
	void processStreamMT(){ // this should be MT，用于处理输入流
		int size;//从文件中读取的数据大小
		int sizeprev;//上一次数据大小
		int bufidx =0;//缓冲区索引，
		int hashbufidx = 0;//哈希缓冲区索引，用于交替使用两个哈希缓冲区
		int start=0;
		bool flag = true;
		ulong_t progress = 0;
		const ulong_t prog_delta = (ulong_t) (FileRead::fsize / 20);//进度更新间隔
		ulong_t nextupd = prog_delta;//下一次进度更新位置
		int progcntr = 0; //进度计数器
          
          size = (unsigned int) loadStream(seqbufs[bufidx]);//从文件中加载数据,offset
		bufAdjust(seqbufs[bufidx], start, size-1);//
		//此时的状态idx=A，readBeg=idx
		streamHashPart(seqbufs[bufidx],start, size - 1, hashbufidx);//存储k-mer的函数
		
		sizeprev = filep.getReadsize();
		start = 0;
		bufidx = (bufidx + 1) & 1;//切换到下一个缓冲区
		size = (unsigned int) loadStream(seqbufs[bufidx]);//读取文件到另一个缓冲区
		while(flag){
			if(size == 0) {
				CS.storeValsInSketchMT(hashbufidx);
				flag = false;
				cout << "\r100%\n" << flush;
				break;
			}

			
			thrd_data[0].id = 0;
			thrd_data[0].hashbufidx = hashbufidx;
			thrd_data[0].cs = &CS;
			pthread_create(&thrds[0], NULL, storeValsWrap, &thrd_data[0]);//哈希值求出来了，存到sktech中了


			hashbufidx  = (hashbufidx + 1) & 1;//toggle
			KmerHash::resetKMBuf(hashbufidx);
			thrd_data[1].id = 1;
			thrd_data[1].buf = seqbufs[bufidx];
			thrd_data[1].start = start;
			thrd_data[1].end = size-1;
			thrd_data[1].hashbufidx = hashbufidx;
			pthread_create(&thrds[1], NULL, streamHashPartWrap, &thrd_data[1]);//保存k-mer

			bufidx = (bufidx + 1) & 1;//toggle
			thrd_data[2].buf = seqbufs[bufidx];
			thrd_data[2].id = 2;

			pthread_create(&thrds[2], NULL, loadStreamWrap, &thrd_data[2]); //从文件中加载数据

			void * status;
			pthread_join(thrds[0], &status);

			pthread_join(thrds[1], &status);

			pthread_join(thrds[2], &status);
			progress += sizeprev;
			if(progress >= nextupd){
				progcntr += 5;
				cout << "\r" << progcntr << "%" << flush;
				nextupd += prog_delta;
			}
			sizeprev = filep.getReadsize();
			size = thrd_data[2].size;

		}
	}



	void processStream(){ // this should be MT
		 int size;
		 int start=0;


		while((size = (  int) loadStream(seqbufs[0])) > 0){
			start = 0;
			bufAdjust(seqbufs[0], start, size-1);
			streamHashPart(seqbufs[0],start, size - 1, 0);


			CS.storeValsInSketch(0);

			KmerHash::resetKMBuf(0);
		}
	}

	static inline void streamHashPart(byte_t * & seqbuf, int  start, int  end, int hashbufidx){//hash a large part of the stream - MT 处理序列数据的一部分，用于哈希化
		int ret= _valX;//15
		int idx = start;
		partAdjust(seqbuf, start, end, idx);//这不是空函数吗

		while(ret != _valEOFBUF){//9
			ret = KmerHash::lShiftSeqHashPart(seqbuf, start, end, idx, hashbufidx );//就是将k-mer存起来

			if(ret == _valEOFBUF) break;
			ret = streamSkip(seqbuf, start, end, idx);
			if(idx > end) break;
		}
		
	}

	static INLINE void skipCommentFasta(byte_t * & seqbuf,  int &idx,  int &end){
		int val;
		bool sawNL = (KmerHash::last_segment_status == ON_COMMENT_LINE_NL) ? true : false;
		while(idx <= end){
			val = Kmer::bpval(seqbuf[idx]);
			if((val == _valNL) || (val == _valCR)){
				sawNL = true;
				KmerHash::last_segment_status = ON_COMMENT_LINE_NL;
			}
			if(sawNL && ((val != _valNL) && (val != _valCR))) {
				KmerHash::last_segment_status = SEG_VALID;
				break;
			}
			++idx;
		}
	}

	static INLINE void bufAdjustFasta (byte_t * & seqbuf,  int & start,  int  end){
		if(KmerHash::last_segment_status == SEG_VALID) return;
		skipCommentFasta(seqbuf, start, end);
	}
	static inline void partAdjustFasta (byte_t * & seqbuf, int & start, int & end, int & idx ){}



	static inline void bufAdjustFastQ (byte_t * & seqbuf,  int & start,  int  end ){
		if((KmerHash::last_segment_status == SEG_VALID) || (KmerHash::last_segment_status == SEEN_VALID_NL)) {//上一个代表reads后的换行或者标识符
			readBeg -= Constants::seq_buf_size; //从0开始
			if( KmerHash::last_segment_status == SEG_VALID) return;
		}
		skipFastQLines(seqbuf, start, end);
	}

	static inline void partAdjustFastQ (byte_t * & seqbuf, int & start, int & end, int & idx ){}

	static INLINE int streamSkipFasta (byte_t * & seqbuf, int & start, int & end, int & idx){ 
		int val = Kmer::bpval(seqbuf[idx]);

		if((val == _valNL)|| (val == _valCR)) {
			++idx;
			if(idx <= end){
				val = Kmer::bpval(seqbuf[idx]);

				if((idx <= end) && ((val == _valNL)|| (val == _valCR))) {
					++idx;
				}
			}

			return 0;
		}
		if((val == _valN) || (val== _valX)){//reset on symbols such as N, U, M, S etc..

			KmerHash::restartKMbufNewSeg();
			++idx;
		}else if(val == _valGT){
			KmerHash::restartKMbufNewSeg();
			KmerHash::last_segment_status = ON_COMMENT_LINE;
			++idx;
			skipCommentFasta(seqbuf, idx, end);
		}

		return 1;
	}

	static inline void skipFastQLines(byte_t * & seqbuf, int & idx, int & end){//跳过多余的行
		int val;

		if((KmerHash::last_segment_status == SEEN_VALID_NL) && ((idx + idx - readBeg + 2) <= end )){ 
			idx += (idx - readBeg + 2);
			KmerHash::last_segment_status = ON_COMMENT_LINE;
		}else{
			if(idx > end) return;

			int val = Kmer::bpval(seqbuf[idx]);


			if(KmerHash::last_segment_status == SEEN_VALID_NL){
				KmerHash::last_segment_status = FQ_SEEN_PLUS;
				++idx;
				if(idx > end) return;
			}


			if(KmerHash::last_segment_status == FQ_SEEN_PLUS){

				KmerHash::last_segment_status = FQ_SEEN_PLUS_NL;
				++idx;
				if(idx > end) return;
			}


			if(KmerHash::last_segment_status == FQ_SEEN_PLUS_NL){

				KmerHash::last_segment_status = FQ_IN_QUAL;
			}


			if(KmerHash::last_segment_status == FQ_IN_QUAL){
				while(idx <= end){
					val = Kmer::bpval(seqbuf[idx]);
					if(val == _valNL){
						KmerHash::last_segment_status = FQ_SEEN_QUAL_NL;
						++idx;
						break;
					}
					++idx;
				}
				if(idx > end) return;
			}


			if(KmerHash::last_segment_status == FQ_SEEN_QUAL_NL){//get to next valid
				KmerHash::last_segment_status = ON_COMMENT_LINE;
			}



		}

		if(KmerHash::last_segment_status == ON_COMMENT_LINE){//开头部分，
			while(idx <= end){//处理标识符
				val = Kmer::bpval(seqbuf[idx]);
				if(val == _valNL){//换行
					KmerHash::last_segment_status = SEG_VALID;//标识符
					++idx;
					readBeg = idx;
					break;
				}
				++idx;
			}
			if(idx > end) return;
		}

	}


	static inline int streamSkipFastaQ (byte_t * & seqbuf, int & start, int & end, int & idx){ 
		int val = Kmer::bpval(seqbuf[idx]);
		if((val == _valNL)) {
			KmerHash::last_segment_status = SEEN_VALID_NL;
			KmerHash::restartKMbufNewSeg();
			++idx;
			skipFastQLines(seqbuf, idx, end);
			return 0;
		}
		KmerHash::restartKMbufNewSeg();
		++idx;
		return 1;
	}

	ulong_t * processAndEstimate(int num_freq){//处理数据并估计k-mer丰度，返回一个数组，包含估计的F0和fi值

		do{
			cout << "Processing " << filep.getfname() << " ...\n";
			if(Constants::streamProcessMT){
				processStreamMT();//问题在这
			}else{
				processStream();
			}

			KmerHash::last_segment_status = ON_COMMENT_LINE;
			KmerHash::resetKMBuf(0);
			KmerHash::resetKMBuf(1);
			KmerHash::restartKMbufNewSeg();

		}while(filep.hasMoreFile());//

		CS.setnfreq(num_freq);



		if(Constants::streamProcessMT){
			CS.analyzeSketchMT();//找最佳w，存起来
		}else{
			CS.analyzeSketch();
		}

		return CS.computeAllF(num_freq);//计算根据公式

	}
};


byte_t * KmerLight::seqbufs[nbufs];
Rands KmerLight::rands[Constants::copies];

void (*KmerLight::partAdjust) (byte_t * & seqbuf, int & start, int & end, int & idx );
void (*KmerLight::bufAdjust) (byte_t * & seqbuf,  int & start,  int  end);

int (*KmerLight::streamSkip) (byte_t * & seqbuf, int & start, int & end, int & idx);

pthread_t KmerLight::thrds[3];
KMLThrdData KmerLight::thrd_data[3];
unsigned long KmerLight::memsize=0;
int KmerLight::readBeg = -1;


FileRead KmerLight::filep;



#endif /* SRC_KMERLIGHT_H_ */

using System;
using System.Numerics;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        //Stopwatch stftTimer = new Stopwatch();
        //public float[][] timeFreqData;
        //public int wSamp;
        //public static double stftWatch;
        //public Complex[] twiddles;
        public float[][] timeFreqData;
        public int wSamp;
        private float[][] Y;
        private int N;
        private Complex[] xx;
        //private FFT fft = new FFT();
        private Complex[] twiddles;
        private float fftMax;
        private int size;
        private int blockSize;
        private Stopwatch sw = new Stopwatch();
        public static double timeTakenSTFT;

        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];

            // Loop 1 storing data in twiddles
            for (ii = 0; ii < wSamp; ii++)
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            }

            timeFreqData = new float[wSamp/2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            // loop 2 Twiddles are not used in this loop
            // data stored in compX 
            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest /wSamp;
            // loop 3 No compX is used in this loop
            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            timeFreqData = stft(compX, wSamp);

            	
        } // All there is no data dependencies between three loops

        //float[][] stft(Complex[] x, int wSamp)
        //{
        //    stftTimer.Start();

        //    //int ii = 0;
        //    //int jj = 0;
        //    int kk = 0;
        //    int ll = 0;
        //    int N = x.Length;
        //    float fftMax = 0;

        //    float[][] Y = new float[wSamp / 2][];

        //    for (ll = 0; ll < wSamp / 2; ll++)
        //    {
        //        Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
        //    }


        //    // Parallel.For(0, 2* Math.Floor((double)N/(double)wSamp) - 1 , () => 0f, 

        //    //Thread implementation

        //   Thread[] fftThreads = new Thread[MainWindow.numThreads];
        //   for (int i = 0; i<MainWindow.numThreads; i++)

        //    for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
        //    {
        //        Complex[] temp = new Complex[wSamp];
        //        Complex[] tempFFT = new Complex[wSamp];

        //        for (int jj = 0; jj < wSamp; jj++)
        //        {
        //            temp[jj] = x[ii * (wSamp / 2) + jj];
        //        }


        //        tempFFT = fft(temp);

        //        for (kk = 0; kk < wSamp / 2; kk++)
        //        {
        //            Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

        //            if (Y[kk][ii] > fftMax)
        //            {
        //                fftMax = Y[kk][ii];
        //            }
        //        }


        //    }

        //    for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
        //    {
        //        for (kk = 0; kk < wSamp / 2; kk++)
        //        {
        //            Y[kk][ii] /= fftMax;
        //        }
        //    }

        //    stftTimer.Stop();
        //    stftWatch = stftTimer.Elapsed.TotalSeconds;
        //    return Y;
        //}

        float[][] stft(Complex[] x, int wSamp)
        {

            N = x.Length;
            xx = x;
            fftMax = 0;
            size = 2 * (int)Math.Floor(N / (double)wSamp);
            blockSize = (size - 1 + MainWindow.numThreads - 1) / MainWindow.numThreads;

            Y = new float[wSamp / 2][];

            Parallel.For(0, wSamp / 2, new ParallelOptions { MaxDegreeOfParallelism = MainWindow.numThreads }, ll =>
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            });

           

            Thread[] fftThreads = new Thread[MainWindow.numThreads];

            for (int i = 0; i < MainWindow.numThreads; i++)
            {
                int id = i;
                fftThreads[i] = new Thread(freqSTFT);
                fftThreads[i].Start(i);
            }
            for (int j = 0; j < MainWindow.numThreads; j++)
            {
                fftThreads[j].Join();
            }


            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
        }

        public void freqSTFT(object threadID)
        {
            int id = (int)threadID;
            int start = id * blockSize;
            int end = Math.Min(start + blockSize, size - 1);
            // dependencies
            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            for (int ii = start; ii < end; ii++)
            {
                for (int jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = xx[ii * (wSamp / 2) + jj];
                }

                //tempFFT = FFT.IterativeFFT(temp, wSamp, twiddles);
                tempFFT = fft(temp);

                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }
            }
        }

        Complex[] fft(Complex[] x)
        {
            int ii = 0;
            int kk = 0;
            int N = x.Length;

            Complex[] Y = new Complex[N];

            // NEED TO MEMSET TO ZERO?

            if (N == 1)
            {
                Y[0] = x[0];
            }
            else{

                Complex[] E = new Complex[N/2];
                Complex[] O = new Complex[N/2];
                Complex[] even = new Complex[N/2];
                Complex[] odd = new Complex[N/2];

                for (ii = 0; ii < N; ii++)
                {

                    if (ii % 2 == 0)
                    {
                        even[ii / 2] = x[ii];
                    }
                    if (ii % 2 == 1)
                    {
                        odd[(ii - 1) / 2] = x[ii];
                    }
                }

                E = fft(even);
                O = fft(odd);

                for (kk = 0; kk < N; kk++)
                {
                    Y[kk] = E[(kk % (N / 2))] + O[(kk % (N / 2))] * twiddles[kk * wSamp / N];
                }
            }

           return Y;
           
        }

        
    }
}

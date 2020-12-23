using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Threading;

namespace Graphic
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {

        DispatcherTimer timer = new DispatcherTimer();

        public MainWindow()
        {
            InitializeComponent();
        }

        private Model3DGroup MainModel3Dgroup = new Model3DGroup();

        // The camera.
        private PerspectiveCamera _theCamera;
        private double a;
        private double b;
        private double c;
        private double percent;
        private int count = 30;
        private int dCount;
        private double u;
        private bool go;

        
        const double dt = 0.02;

        // Make the surface's points and triangles.
        const double xMin = -2;
        const double xMax = 2;
        const double dx = 0.5;
        const double zMin = -2;
        const double zMax = 2;
        const double dz = 0.5;

        private double[] x, z;
        Complex[][] psi; // Волновая функция
        Complex[][] alfaX; // Коэффициенты для метода сетчатой прогонки
        Complex[][] betaX;
        Complex[][] alfaZ;
        Complex[][] betaZ;

        Complex CoefAx()
        {
            return new Complex(0, -(dt / (4 * dx * dx)));
        }

        Complex CoefAz()
        {
            return new Complex(0, -(dt / (4 * dz * dz)));
        }

        //Функция для расчета коэффициента B
        Complex CoefBx()
        {
            return new Complex(0, -(dt / (4 * dx * dx)));
        }

        Complex CoefBz()
        {
            return new Complex(0, -(dt / (4 * dz * dz)));
        }

        //Функция для расчета коэффициента C
        Complex CoefCx(int i, int j)
        {
            return new Complex(1, dt / (2*dx * dx) + U(i, j) * dt / 4);
        }

        Complex CoefCz(int i, int j)
        {
            return new Complex(1, dt / (2*dz * dz) + U(i, j) * dt / 4);
        }

        //Функция для расчета коэффициента D
        Complex CoefDx(int i, int j, Complex[][] psi)
        {
            //return new Complex(
            //    psi[i][j].Real + psi[i][j].Imaginary * dt * (U(i, j) / 2 + 2 / (dz * dz) + 1 / (dx * dx)) - dt / (2 * dx * dx) * (psi[i + 1][j].Imaginary + psi[i - 1][j].Imaginary) -
            //    dt / (dz * dz) * (psi[i][j + 1].Imaginary + psi[i][j - 1].Imaginary),
            //    psi[i][j].Imaginary + -psi[i][j].Real * dt * (U(i, j) / 2 + 2 / (dz * dz) + 1 / (dx * dx)) + dt / (2 * dx * dx) * (psi[i + 1][j].Real + psi[i - 1][j].Real) +
            //    dt / (dz * dz) * (psi[i][j + 1].Real + psi[i][j - 1].Real));
            return psi[i][j]*Complex.One - psi[i][j]*Complex.ImaginaryOne*dt*(U(i,j)/4 + 1/(dz*dz) + 1/(2*dx*dx)) + Complex.ImaginaryOne*dt/(4*dx*dx)*(psi[i+1][j]+ psi[i - 1][j])+ Complex.ImaginaryOne * dt / (2 * dz * dz) * (psi[i][j+1] + psi[i][j-1]);
        }

        Complex CoefDz(int i, int j, Complex[][] psi)
        {
            //return new Complex(
            //    psi[i][j].Real + psi[i][j].Imaginary * dt * (U(i, j) / 2 + 2 / (dx * dx) + 1 / (dz * dz)) - dt / (2 * dz * dz) * (psi[i][j + 1].Imaginary + psi[i][j - 1].Imaginary) -
            //    dt / (dx * dx) * (psi[i + 1][j].Imaginary + psi[i - 1][j].Imaginary),
            //    psi[i][j].Imaginary - psi[i][j].Real * dt * (U(i, j) / 2 + 2 / (dx * dx) + 1 / (dz * dz)) + dt / (2 * dz * dz) * (psi[i][j + 1].Real + psi[i][j - 1].Real) +
            //    dt / (dx * dx) * (psi[i + 1][j].Real + psi[i - 1][j].Real));
            return psi[i][j] * Complex.One - psi[i][j] * Complex.ImaginaryOne * dt * (U(i, j) / 4 + 1 / (dx * dx) + 1 / (2 * dz * dz)) + Complex.ImaginaryOne * dt / (2 * dx * dx) * (psi[i + 1][j] + psi[i - 1][j]) + Complex.ImaginaryOne * dt / (4 * dz * dz) * (psi[i][j + 1] + psi[i][j - 1]);
        }

        // Рассчёт коэффициента альфа
        Complex FuncAlfaX(int i, int j, Complex prevAlfa)
        {
            return -CoefBx() / (CoefCx(i, j) + CoefAx() * prevAlfa);
        }

        Complex FuncAlfaZ(int i, int j, Complex prevAlfa)
        {
            return -CoefBz() / (CoefCz(i, j) + CoefAz() * prevAlfa);
        }



        // Рассчёт коффициента бетта
        Complex FuncBetaX(int i, int j, Complex[][] psi, Complex prevAlfa, Complex prevBeta)
        {
            return (CoefDx(i, j, psi) - CoefAx() * prevBeta) / (CoefCx(i, j) + CoefAx() * prevAlfa);
        }

        Complex FuncBetaZ(int i, int j, Complex[][] psi, Complex prevAlfa, Complex prevBeta)
        {
            return (CoefDz(i, j, psi) - CoefAz() * prevBeta) / (CoefCz(i, j) + CoefAz() * prevAlfa);
        }

        double U(int i, int j)
        {
            return (Math.Abs(x[i]) > a || Math.Abs(z[j]) > b) ? 0 : -u;
        }

        // The camera's current location.
        private double _cameraPhi = Math.PI / 6.0; // 30 degrees
        private double _cameraTheta = Math.PI / 6.0; // 30 degrees

        private double _cameraR = 3.0;
        //private double CameraR = 13.0;

        // The change in CameraPhi when you press the up and down arrows.
        private const double CameraDPhi = 0.1;

        // The change in CameraTheta when you press the left and right arrows.
        private const double CameraDTheta = 0.1;

        // The change in CameraR when you press + or -.
        private const double CameraDr = 0.1;

        // Create the scene.
        // MainViewport is the Viewport3D defined
        // in the XAML code that displays everything.
        private void WindowLoaded(object sender, RoutedEventArgs e)
        {
            // Give the camera its initial position.
            _theCamera = new PerspectiveCamera {FieldOfView = 120};
            MainViewport.Camera = _theCamera;
            PositionCamera();
            
            u = 3;
        }

        void Init() // Функция инцилизации компонентов
        {
            var xCount = (int) (xMax * 2 / dx) + 1;
            var zCount = (int) (zMax * 2 / dz) + 1;
            psi = new Complex[xCount][];

                for (int j = 0; j < xCount; j++)
                {
                    psi[j] = new Complex[zCount];
                }

            x = new double[xCount];
            z = new double[zCount];

            dCount = 0;
            alfaX = new Complex[xCount][];
            betaX = new Complex[xCount][];
            alfaZ = new Complex[xCount][];
            betaZ = new Complex[xCount][];
            for (int i = 0; i < xCount; i++)
            {
                alfaX[i] = new Complex[zCount];
                betaX[i] = new Complex[zCount];
                alfaZ[i] = new Complex[zCount];
                betaZ[i] = new Complex[zCount];
            }

            for (int i = 0; i <= (xMax - xMin) / dx; ++i)
                x[i] = i * dx + xMin;
            for (int i = 0; i <= (zMax - zMin) / dz; ++i)
                z[i] = i * dz + zMin;

        }


        // Define the lights.
        private void DefineLights()
        {
            var ambientLight = new AmbientLight(Colors.Gray);
            var directionalLight = new DirectionalLight(Colors.Gray, new Vector3D(-1.0, -3.0, -2.0));
            MainModel3Dgroup.Children.Add(ambientLight);
            MainModel3Dgroup.Children.Add(directionalLight);
        }

        // Add the model to the Model3DGroup.
        private void DefineModel(Model3DGroup modelGroup, int t)
        {

            // Define lights.
            DefineLights();
            // Make a mesh to hold the surface.
            var mesh = new MeshGeometry3D();

            for (int i = 0; i < x.Length-1; ++i)
            {
                for (int j = 0; j < z.Length - 1; ++j)
                {
                    // Make points at the corners of the surface
                    // over (x, z) - (x + dx, z + dz).

                    var p00 = new Point3D(x[i], psi[i][j].Magnitude, z[j]);
                    var p10 = new Point3D(x[i + 1], psi[i + 1][j].Magnitude, z[j]);
                    var p01 = new Point3D(x[i], psi[i][j + 1].Magnitude, z[j + 1]);
                    var p11 = new Point3D(x[i + 1], psi[i + 1][j + 1].Magnitude, z[j + 1]);

                    // Add the triangles.
                    AddTriangle(mesh, p00, p01, p11);
                    AddTriangle(mesh, p00, p11, p10);
                }
            }

            Console.WriteLine(mesh.Positions.Count + " points");
            Console.WriteLine(mesh.TriangleIndices.Count / 3 + " triangles");

            // Make the surface's material using a solid orange brush.
            var surfaceMaterial = new DiffuseMaterial(Brushes.Orange);

            // Make the mesh's model. Make the surface visible from both sides.
            var surfaceModel = new GeometryModel3D(mesh, surfaceMaterial) {BackMaterial = surfaceMaterial};
            
            // Add the model to the model groups.
            modelGroup.Children.Add(surfaceModel);
        }

        // The function that defines the surface we are drawing.
        private double F(int i, int j)
        {
            var f = 0.5 * Math.Exp(-(x[i] * x[i]) / (2 * c * c)) * Math.Exp(-(z[j] * z[j]) / (2 * c * c));
            return f;
        }

        // Add a triangle to the indicated mesh.
        private void AddTriangle(MeshGeometry3D mesh, Point3D point1, Point3D point2, Point3D point3)
        {
            // Get the points' indices.
            var index1 = AddPoint(mesh.Positions, point1);
            var index2 = AddPoint(mesh.Positions, point2);
            var index3 = AddPoint(mesh.Positions, point3);

            // Create the triangle.
            mesh.TriangleIndices.Add(index1);
            mesh.TriangleIndices.Add(index2);
            mesh.TriangleIndices.Add(index3);
        }

        // Create the point and return its new index.
        private static int AddPoint(Point3DCollection points, Point3D point)
        {
            // Create the point and return its index.
            points.Add(point);
            return points.Count - 1;
        }

        // Adjust the camera's position.
        private void WindowKeyDown(object sender, KeyEventArgs e)
        {
            switch (e.Key)
            {
                case Key.Up:
                    _cameraPhi += CameraDPhi;
                    if (_cameraPhi > Math.PI / 2.0) _cameraPhi = Math.PI / 2.0;
                    break;
                case Key.Down:
                    _cameraPhi -= CameraDPhi;
                    if (_cameraPhi < -Math.PI / 2.0) _cameraPhi = -Math.PI / 2.0;
                    break;
                case Key.Left:
                    _cameraTheta += CameraDTheta;
                    break;
                case Key.Right:
                    _cameraTheta -= CameraDTheta;
                    break;
                case Key.Add:
                case Key.OemPlus:
                    _cameraR -= CameraDr;
                    if (_cameraR < CameraDr) _cameraR = CameraDr;
                    break;
                case Key.Subtract:
                case Key.OemMinus:
                    _cameraR += CameraDr;
                    break;
            }

            // Update the camera's position.
            PositionCamera();
        }

        // Position the camera.
        private void PositionCamera()
        {
            // Calculate the camera's position in Cartesian coordinates.
            double y = _cameraR * Math.Sin(_cameraPhi);
            double hyp = _cameraR * Math.Cos(_cameraPhi);
            double x = hyp * Math.Cos(_cameraTheta);
            double z = hyp * Math.Sin(_cameraTheta);
            _theCamera.Position = new Point3D(x, y, z);

            // Look toward the origin.
            _theCamera.LookDirection = new Vector3D(-x, -y, -z);

            // Set the Up direction.
            _theCamera.UpDirection = new Vector3D(0, 1, 0);

            // Console.WriteLine("Camera.Position: (" + x + ", " + y + ", " + z + ")");
        }

        private void RunButtonClick(object sender, RoutedEventArgs e)
        {
            Init();
            var aCond = double.TryParse(ATextBox.Text, out a);
            var bCond = double.TryParse(BTextBox.Text, out b);
            var cCond = double.TryParse(CTextBox.Text, out c);

            if (!aCond && !bCond && !cCond)
            {
                MessageBox.Show("Проверьте введённые значения");
                return;
            }

            go = true;
            // Add the group of models to a ModelVisual3D.
            var modelVisual = new ModelVisual3D { Content = MainModel3Dgroup };

            // Add the main visual to the viewportt.
            MainViewport.Children.Add(modelVisual);
            ATextBox.IsEnabled = false;
            BTextBox.IsEnabled = false;
            CTextBox.IsEnabled = false;
            timer.Tick += new EventHandler(timer_Tick);
            timer.Interval = new TimeSpan(0, 0, 1);
            timer.Start();
        }

        private void StopButtonClick(object sender, RoutedEventArgs e)
        {
            ATextBox.IsEnabled = true;
            BTextBox.IsEnabled = true;
            CTextBox.IsEnabled = true;

            go = false;
            MainModel3Dgroup.Children.Clear();
            MainViewport.Children.Clear();
        }


        private void EvolutionButtonClick(object sender, RoutedEventArgs e)
        {
            timer.Tick += new EventHandler(timer_Tick);
            timer.Interval = TimeSpan.FromMilliseconds(10);
            timer.Start();
           
        }


        void Generate()
        {
            if (dCount == 0)
                for (int i = 0; i < x.Length; ++i)
                {
                    for (int j = 0; j < z.Length; ++j)
                    {
                        psi[i][j] = F(i, j);
                    }
                }
            else
            {
                for (int i = 0; i < x.Length - 1; ++i)
                {
                    //начальные условия
                    for (int j = 0; j < z.Length-1; ++j)
                    {
                        // Прямой ход сетчатой прогонки
                        if (i == 0 || j == 0 || i == 1 || j == 1)
                        {
                            alfaX[i][j] = 0;
                            betaX[i][j] = 0;
                        }
                        else
                        {
                            alfaX[i][j] = (FuncAlfaX(i, j - 1, alfaX[i ][j - 1]));
                            betaX[i][j] = (FuncBetaX(i, j-1, psi, alfaX[i ][j - 1], betaX[i ][j - 1]));
                        }
                    }

                    for (int j = z.Length - 1; j > -1; --j)
                    {
                        if (i == 0 || j == 0 || i == x.Length - 1 || j == z.Length - 1)
                        {
                            psi[i][j] = 0;
                        }
                        else
                        {
                            psi[i][j] = alfaX[i][j + 1] * psi[i][j + 1] + betaX[i][j + 1];
                        }
                    }
                }

                for (int i = 0; i < x.Length - 1; ++i)
                {
                    for (int j = 0; j < z.Length - 1; ++j)
                    {
                        // Прямой ход сетчатой прогонки
                        if (i == 0 || j == 0 || i ==1 || j ==  1)
                        {
                            alfaZ[i][j] = 0;
                            betaZ[i][j] = 0;
                        }
                        else
                        {
                            alfaZ[j][i] = (FuncAlfaZ(j-1, i, alfaZ[j - 1][i]));
                            betaZ[j][i] = (FuncBetaZ(j-1, i, psi, alfaZ[j - 1][i], betaZ[j - 1][i]));
                        }

                    }

                    for (int j = z.Length - 1; j > -1; --j)
                    {
                        if (i == 0 || j == 0 || i == x.Length - 1 || j == z.Length - 1)
                        {
                            psi[i][j] = 0;
                        }
                        else
                        {
                            psi[j][i] = alfaZ[j + 1][i] * psi[j+1][i] + betaZ[j + 1][i];
                        }
                    }
                }
            }
        }

        private void timer_Tick(object sender, EventArgs e)
        {
            // Create the model.

            if (dCount > 0)
            {
                MainModel3Dgroup.Children.Clear();
                //MainViewport.Children.Clear();
            }
            Generate();
            DefineModel(MainModel3Dgroup, dCount);
            

            ++dCount;
                if (!go)
                timer.Stop();
        }

    }
}
        




using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
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

        // The camera's current location.
        private double _cameraPhi = Math.PI / 6.0;       // 30 degrees
        private double _cameraTheta = Math.PI / 6.0;     // 30 degrees

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
        private void DefineModel(Model3DGroup modelGroup)
        {
            // Define lights.
            DefineLights();
            // Make a mesh to hold the surface.
            var mesh = new MeshGeometry3D();

            // Make the surface's points and triangles.
            const double xMin = -1.5;
            const double xMax = 1.5;
            const double dx = 0.05;
            const double zMin = -1.5;
            const double zMax = 1.5;
            const double dz = 0.05;

            //const double xmin = -5;
            //const double xmax = 5;
            //const double dx = 0.5;
            //const double zmin = -5;
            //const double zmax = 5;
            //const double dz = 0.5;

            for (double x = xMin; x <= xMax - dx; x += dx)
            {
                for (double z = zMin; z <= zMax - dz; z += dx)
                {
                    // Make points at the corners of the surface
                    // over (x, z) - (x + dx, z + dz).
                    var p00 = new Point3D(x, F(x, z), z);
                    var p10 = new Point3D(x + dx, F(x + dx, z), z);
                    var p01 = new Point3D(x, F(x, z + dz), z + dz);
                    var p11 = new Point3D(x + dx, F(x + dx, z + dz), z + dz);

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
        private double F(double x, double z)
        {
            var f = 1.5 * Math.Exp(-(x * x) / (2 * c * c)) * Math.Exp(-(z * z) / (2 * c * c));
            return f;

            //double r2 = x * x + z * z;
            //return 8 * Math.Cos(r2 / 2) / (2 + r2);
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
            var aCond = double.TryParse(ATextBox.Text, out a);
            var bCond = double.TryParse(BTextBox.Text, out b);
            var cCond = double.TryParse(CTextBox.Text, out c);

            if (!aCond && !bCond && !cCond)
            {
                MessageBox.Show("Проверьте введённые значения");
                return;
            }


            ATextBox.IsEnabled = false;
            BTextBox.IsEnabled = false;
            CTextBox.IsEnabled = false;
            // Create the model.
            DefineModel(MainModel3Dgroup);

            // Add the group of models to a ModelVisual3D.
            var modelVisual = new ModelVisual3D { Content = MainModel3Dgroup };

            // Add the main visual to the viewportt.
            MainViewport.Children.Add(modelVisual);
        }

        private void StopButtonClick(object sender, RoutedEventArgs e)
        {
            ATextBox.IsEnabled = true;
            BTextBox.IsEnabled = true;
            CTextBox.IsEnabled = true;

            MainModel3Dgroup.Children.Clear();
            MainViewport.Children.Clear();
        }
        public void Eee()
        {
            timer.Tick += new EventHandler(timer_Tick);
            timer.Interval = new TimeSpan(0, 0, 1);
            timer.Start();
        }

        private void timer_Tick(object sender, EventArgs e)
        {
            // код здесь
        }

    }
}
        




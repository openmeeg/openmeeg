class Openmeeg < Formula
  desc "Package for low-frequency bio-electromagnetism forward models"
  homepage "https://openmeeg.github.io"
  head "https://github.com/openmeeg/openmeeg.git"
  url "https://github.com/openmeeg/openmeeg/archive/2.4.1.tar.gz"
  version "2.4.1"

  devel do
    url "https://github.com/openmeeg/openmeeg/archive/master.zip"
    version "development"
  end

  depends_on "cmake" => :build
  depends_on "hdf5"
  depends_on "libmatio"
  depends_on "zlib" unless OS.mac?
  depends_on "openblas" unless OS.mac?

  needs :openmp

  def install
    args = std_cmake_args + %w[
      -DBLA_VENDOR=OpenBLAS
      -DBUILD_SHARED_LIBS=ON
      -DENABLE_PYTHON=OFF
      -DBUILD_DOCUMENTATION=OFF
      -DUSE_PROGRESSBAR=ON
      -DUSE_OMP=ON
    ]

    mkdir "build" do
      args << ".."
      system "cmake", *args
      system "make", "install"
    end
  end

  test do
    system "#{bin}/om_assemble"
  end
end
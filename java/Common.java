import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Common {

  public static double[] toDouble(String[] rs) throws NumberFormatException {
    double[] d = new double[rs.length];
    int pos = 0;
    for (String s : rs) {
      s = s.trim();
      if (isDouble(s)) {
        d[pos] = Double.parseDouble(s.trim());
        pos++;
      } else {
        throw new NumberFormatException();
      }
    }
    return d;
  }

  public static boolean isDouble(String string) {
    try {
      Double.valueOf(string);
      return true;
    } catch (NumberFormatException e) {
      return false;
    }
  }

  public static double setDecimalState(double d, int n) {
    int x = 1;
    while (n > 0) {
      x *= 10;
      n--;
    }
    return ((int) (d * x)) / (double) x;
  }

  public static String toString(String[] d, String sep) {
    StringBuffer b = new StringBuffer();
    for (String s : d) {
      b.append(s + sep);
    }
    return b.substring(0, b.length() - sep.length());
  }

  public static String toString(double[] d2, String sep) {
    StringBuffer b = new StringBuffer();
    int i = 0;
    for (; i < d2.length - 1; i++) {
      b.append(setDecimalState(d2[i], 2) + sep);
    }
    b.append(setDecimalState(d2[i], 2));
    return b.toString();
  }

  public static BufferedReader getReader(String cbsin)
      throws FileNotFoundException {
    return new BufferedReader(new FileReader(new File(cbsin)));
  }

  public static BufferedWriter getWriter(String out) throws IOException {
    return new BufferedWriter(new FileWriter(new File(out)));
  }

  public static int firstIndexOf(Object o, Object[] array) {
    for (int i = 0; i < array.length; i++) {
      if (array[i].equals(o)) {
        return i;
      }
    }
    return -1;
  }

  public static int argmin(double[] is) {
    double min = Double.POSITIVE_INFINITY;
    int pos = -1;
    for (int idx = 0; idx < is.length; idx++) {
      double i = is[idx];
      if (min > i) {
        min = i;
        pos = idx;
      }
    }
    return pos;
  }
}
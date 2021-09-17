function d = point_on_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm( dot(-a, b) * a );
end
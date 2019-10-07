#include "MyArea.hpp"

bool MyArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr)
{
  Gtk::Allocation allocation = get_allocation();
  const int width = allocation.get_width();
  const int height = allocation.get_height();

  // choose which coordinates to use for 2d image
  const int x_coord = 0;
  const int y_coord = 1;

  // find min and max on each axis
  double x_min = DBL_MAX;
  double x_max = DBL_MIN;
  double y_min = DBL_MAX;
  double y_max = DBL_MIN;
  for (int i = 0; i < c->get_n(); i++)
  {
    if (c->get_point(i).coords[x_coord] < x_min)
      x_min = c->get_point(i).coords[x_coord];

    if (c->get_point(i).coords[x_coord] > x_max)
      x_max = c->get_point(i).coords[x_coord];

    if (c->get_point(i).coords[y_coord] < y_min)
      y_min = c->get_point(i).coords[y_coord];

    if (c->get_point(i).coords[y_coord] > y_max)
      y_max = c->get_point(i).coords[y_coord];
  }

  // plot all points
  for (int i = 0; i < c->get_n(); i++)
  {
    cr->save(); // save current drawing context (opaque black)

    // choose color depending on label
    srand(time(0));
    int l = c->get_point(i).label;
    if (l > 0)
    {
      cr->arc((c->get_point(i).coords[x_coord] - x_min) * width / (x_max - x_min), (c->get_point(i).coords[y_coord] - y_min) * height / (y_max - y_min), 5.0, 0.0, 2.0 * M_PI); // full circle
      cr->set_source_rgba(l % 2, l % 3 / 2.0, l % 4 / 3.0, 0.6);
    }
    else //We show them as small points
    {
      cr->arc((c->get_point(i).coords[x_coord] - x_min) * width / (x_max - x_min), (c->get_point(i).coords[y_coord] - y_min) * height / (y_max - y_min), 2.0, 0.0, 2.0 * M_PI); // full circle

      cr->set_source_rgba(1, 1, 1, 0);
    }
    /*switch (c->get_point(i).label)
    {
    case 0:
      cr->set_source_rgba(1.0, 0.0, 0.0, 0.6); // red, partially translucent
      break;

    case 1:
      cr->set_source_rgba(0.0, 0.0, 0.8, 0.6); // 0.8 blue, partially translucent
      break;

    case 2:
      cr->set_source_rgba(0.0, 1.0, 0.0, 0.6); // green, partially translucent
      break;

    default:
      cr->set_source_rgba(1.0, 1.0, 0.0, 0.6); // yellow, partially translucent
      break;
    }*/

    cr->fill_preserve();
    cr->restore(); // restore drawing context to opaque black
    cr->stroke();
  }

  return true;
}

#ifndef GRAPHER
#define GRAPHER

#include<SFML/Graphics.hpp>
#include<iostream>
#include<string>
#include<cstring>
#include<sstream>
#include<iomanip>

namespace MyColors {
    const sf::Color Primary = sf::Color::White;
    const sf::Color Secondary = sf::Color::Black;
    const sf::Color point1 = sf::Color::Green;
    const sf::Color point2 = sf::Color::Blue;
}

class MyPoint : public sf::Vector2f {
private:
    sf::Color color;
public:
    MyPoint(float x, float y, sf::Color col = MyColors::Primary) : Vector2<float>(x, y), color(col) {}

    friend class Graph;
};

class Graph {
private:
    sf::RenderWindow window;
    int screenWidth;
    int screenHeight;
    sf::Font font;
    std::vector<MyPoint> points;
    std::vector<MyPoint> valPoints;
    std::vector<std::string> legends;
    std::vector<sf::Color> legendColor;
    sf::VertexArray lines;
    double maxX = 10;
    double maxY = 0.1;
    double minY = 0;

public:
    Graph(int sw = 1450, int sh = 1000)
        : window(sf::VideoMode(sw, sh), "A Graph", sf::Style::Default) {
        window.setFramerateLimit(60);
        screenWidth = sw;
        screenHeight = sh;

        if (!font.loadFromFile("Utils/fonts/arial/ARIAL.TTF")) {
            std::cerr << "Failed to load font!" << std::endl;
            return;
        }
    }

    void drawAxes() {
        double originX = 100;
        double originY = screenHeight - 70;

        // Draw X-axis
        sf::VertexArray xAxis(sf::Lines, 2);
        xAxis[0].position = sf::Vector2f(originX, originY);
        xAxis[1].position = sf::Vector2f(screenWidth - 50, originY);
        xAxis[0].color = MyColors::Primary;
        xAxis[1].color = MyColors::Primary;
        window.draw(xAxis);

        // Draw Y-axis
        sf::VertexArray yAxis(sf::Lines, 2);
        yAxis[0].position = sf::Vector2f(originX, 50);
        yAxis[1].position = sf::Vector2f(originX, originY);
        yAxis[0].color = MyColors::Primary;
        yAxis[1].color = MyColors::Primary;
        window.draw(yAxis);

        int yMarkingInterval = 50;


        std::ostringstream oss1;
        oss1 << std::fixed << std::setprecision(2) << (minY);
        std::string minNumber = oss1.str();


        sf::Text label;
        label.setFont(font);
        label.setString(minNumber);
        label.setCharacterSize(15);
        label.setFillColor(MyColors::Primary);
        label.setPosition(originX - 20, originY - 10);
        window.draw(label);

        for (int i = originY - yMarkingInterval; i > 70; i -= yMarkingInterval) {
            sf::VertexArray marking(sf::Lines, 2);
            marking[0].position = sf::Vector2f(originX - 5, i);
            marking[1].position = sf::Vector2f(originX + 5, i);
            marking[0].color = MyColors::Primary;
            marking[1].color = MyColors::Primary;
            window.draw(marking);

            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << (((double)(originY - i) * ((maxY - minY) / (double)(screenHeight - 140))) + minY);
            std::string numberStr = oss.str();

            sf::Text label;
            label.setFont(font);
            label.setString(numberStr);
            label.setCharacterSize(15);
            label.setFillColor(MyColors::Primary);
            label.setPosition(originX - 30, i - 10);
            window.draw(label);
        }

        int xMarkingInterval = 50;

        for (int i = originX + xMarkingInterval; i < (screenWidth - 50); i += xMarkingInterval) {
            sf::VertexArray marking(sf::Lines, 2);
            marking[0].position = sf::Vector2f(i, originY + 5);
            marking[1].position = sf::Vector2f(i, originY - 5);
            marking[0].color = MyColors::Primary;
            marking[1].color = MyColors::Primary;
            window.draw(marking);

            sf::Text label;
            label.setFont(font);
            label.setString(std::to_string(static_cast<int>((double)(i - originX) * (maxX / (double)(screenWidth - 150)))));
            label.setCharacterSize(15);
            label.setFillColor(MyColors::Primary);
            label.setPosition(i, originY + 10);
            window.draw(label);
        }
    }

    void drawLabels(std::string xLabel = "Epochs", std::string yLabel = "Loss") {
        double originX = 100;
        double originY = screenHeight - 70;

        // Draw X-axis label
        sf::Text xAxisLabel;
        xAxisLabel.setString(xLabel);
        xAxisLabel.setFont(font);
        xAxisLabel.setFillColor(MyColors::Primary);
        xAxisLabel.setPosition(screenWidth / 2 - 50, originY + 25);
        xAxisLabel.setCharacterSize(30);
        window.draw(xAxisLabel);

        // Draw Y-axis label
        sf::Text yAxisLabel(yLabel, font, 30);
        yAxisLabel.setFillColor(MyColors::Primary);
        yAxisLabel.setPosition(15, screenHeight / 2);
        yAxisLabel.setCharacterSize(30);
        yAxisLabel.setRotation(-90.0);
        window.draw(yAxisLabel);
    }

    void pushColoredPoints(const std::vector<double>& newPoints, const sf::Color color) {

        for (int i = 0; i < newPoints.size(); i++) {
            points.push_back(MyPoint(i, newPoints[i], color));
        }

    }

    void pushColoredValPoints(const std::vector<double>& newPoints, const sf::Color color) {

        for (int i = 0; i < newPoints.size(); i++) {
            valPoints.push_back(MyPoint(i, newPoints[i], color));
        }

    }

    void pushPoint(const MyPoint& point) {
        points.push_back(point);
    }

    void pushValPoint(const MyPoint& point) {
        points.push_back(point);
    }

    void clearGraph() {
        points.clear();
        valPoints.clear();
    }

    void plotPoints(const vector<MyPoint>& points) {
        double originX = 100;
        double originY = screenHeight - 70;

        sf::VertexArray lines(sf::LinesStrip, points.size());

        int i = 0;
        for (const auto& point : points) {
            if (point.x > maxX) {
                maxX = point.x + point.x * 0.1;
            }

            if (point.y > maxY) {
                maxY = point.y + point.y * 0.1;
            }

            if (point.y < minY) {
                minY = point.y;
            }

            double x = (point.x * (double)(screenWidth - 150) / maxX) + originX;
            double y = originY - ((point.y - minY) * (double)(screenHeight - 140) / (maxY - minY));

            sf::CircleShape pointShape(2); // Radius of 5
            pointShape.setFillColor(point.color);
            pointShape.setPosition(x - 1, y - 1);
            lines[i].position = sf::Vector2f(x, y);
            lines[i].color = point.color;

            i++;

            window.draw(pointShape);
        }

        window.draw(lines);
    }

    void setLegend(const std::vector<std::string>& lgds, const std::vector<sf::Color>& colors) {
        legends = lgds;
        legendColor = colors;
    }

    void displayLegend() {
        sf::Text legend;

        for (int i = 0; i < legends.size(); i++) {
            sf::CircleShape circle(10);

            circle.setFillColor(legendColor[i]);
            circle.setPosition(screenWidth / (6 - i * 4) - 25, 48);

            legend.setString(legends[i]);
            legend.setFont(font);
            legend.setFillColor(legendColor[i]);
            legend.setCharacterSize(30);
            legend.setPosition(screenWidth / (6 - i * 4), 40);

            window.draw(legend);
            window.draw(circle);
        }
    }

    void start(char x[], char y[]) {
        bool run = true;

        while (run) {
            sf::Event event;
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    run = false;
                }
                else if (event.type == sf::Event::Resized) {
                    sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
                    window.setView(sf::View(visibleArea));
                }
            }

            screenWidth = window.getSize().x;
            screenHeight = window.getSize().y;

            window.clear(MyColors::Secondary);

            drawAxes();
            drawLabels(x, y);
            displayLegend();

            if (points.size() > 0) {
                plotPoints(points);
            }

            if (valPoints.size() > 0) {
                plotPoints(valPoints);
            }

            window.display();
        }
    }
};

#endif

// int main() {
//     char x[] = "Epochs";
//     char y[] = "Loss";

//     // MyPoint point(3, 4, MyColors::point1);

//     std::vector<MyPoint> points = {
//         MyPoint(500, 300, MyColors::point1),
//         {500, 500},
//     };

//     Graph plotter;

//     plotter.setPoints(points);

//     plotter.pushPoint(MyPoint(100, 100, MyColors::point2));

//     plotter.start(x, y);

//     return 0;
// }

// Copyright (C) 2026  Cody Raskin
#pragma once

#include <iostream>

#include <cmath> // for round

#include "colors.hh"

void ProgressBar(double pct, std::string text) {
    // Define the length of the progress bar
    const int barWidth = 20;

    // Calculate the number of filled positions
    int filledLength = static_cast<int>(round(pct * barWidth));

    // Create the progress bar string
    std::string filled, remainder;
    for (int i = 0; i < filledLength; ++i) {
        filled += "=";
    }
    for (int i = filledLength + 1; i < barWidth; ++i) {
        remainder += "·";
    }

    std::string bar = Color::colorize("[", Color::dim());
    bar += Color::colorize(filled, Color::green());
    bar += Color::colorize(">", Color::brightGreen());
    bar += Color::colorize(remainder, Color::dim());
    bar += Color::colorize("]", Color::dim());

    bar += Color::colorize(" " + std::to_string(pct*100) + "%", Color::brightWhite());

    bar += " " + text;

    // Print the progress bar
    std::cout << bar << std::endl;
}
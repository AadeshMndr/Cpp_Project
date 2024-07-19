#include<iostream>
#include<raylib.h>
#include <string>
#include <chrono>
#include <thread>
#include "utils.h"

const int screen_width=1920/2;
const int screen_height=1080/2;


int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(screen_width, screen_height, "Derivative Calculator");
    SetTargetFPS(144);
    // Title
    const std::string title = "Derivative Calculator";
    SetWindowTitle(title.c_str());
    int titleSize = 50;
    // Calculate the position to center the text
    int textWidth = MeasureText(title.c_str(), titleSize);
    int xPosition = (screen_width - textWidth) / 2;
    int yPosition = (screen_height - titleSize) / 2;

    // equation input box
    const int equtionInputBoxWidth = 760;
    const int equtionInputBoxHeight = 40;
    Rectangle equtionInputBox = { (screen_width - equtionInputBoxWidth) / 2, (screen_height - equtionInputBoxHeight) / 3, 760, 40 }; // Define the input box rectangle
    bool equationEditMode = false;  // Flag to indicate if we're editing the input text
    char equationText[64] = "\0";  // Initialize empty input field
    int equationLetterCount = 0;

    // equation input box
    const int variableInputBoxWidth = 760;
    const int variableInputBoxHeight = 40;
    Rectangle variableInputBox = { (screen_width - variableInputBoxWidth) / 2, (screen_height - variableInputBoxHeight) / 1.8, 760, 40 }; // Define the input box rectangle
    bool variableEditMode = false;  // Flag to indicate if we're editing the input text
    char variableText[64] = "\0";  // Initialize empty input field
    int varibaleLetterCount = 0;
    float variableFloat = 0.0f;

    // Backspace key repeat variables
    bool backspaceHeld = false;
    float backspaceTimer = 0.0f;
    const float BACKSPACE_REPEAT_DELAY = 0.5f; // Adjust as needed
    const float BACKSPACE_REPEAT_INTERVAL = 0.05f; // Interval between repeats

    while (!WindowShouldClose())
    {
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            Vector2 mousePos = GetMousePosition();
            // Check if mouse is inside the input box
            if (CheckCollisionPointRec(mousePos, equtionInputBox)) {
                equationEditMode = true; // Enable edit mode
            } else {
                equationEditMode = false; // Disable edit mode if clicked outside
            }
             if (CheckCollisionPointRec(mousePos, variableInputBox)) {
                variableEditMode = true; // Enable edit mode
            } else {
                variableEditMode = false; // Disable edit mode if clicked outside
            }
        }

            if (equationEditMode && equationLetterCount < 64) {
                int key = GetCharPressed();
                // Check for continuous backspace key press
            if (IsKeyDown(KEY_BACKSPACE)) {
                backspaceTimer += GetFrameTime();
                if (!backspaceHeld && backspaceTimer >= BACKSPACE_REPEAT_DELAY) {
                    backspaceHeld = true;
                    backspaceTimer = 0.0f;

                    if (equationLetterCount > 0) {
                        equationLetterCount--;
                        equationText[equationLetterCount] = '\0'; // Clear the character at the current position
                    }
                }
                else if (backspaceHeld && backspaceTimer >= BACKSPACE_REPEAT_INTERVAL) {
                    backspaceTimer = 0.0f;

                    if (equationLetterCount > 0) {
                        equationLetterCount--;
                        equationText[equationLetterCount] = '\0'; // Clear the character at the current position
                    }
                }
            }
            else {
                backspaceHeld = false;
                backspaceTimer = 0.0f;
            }

            if (IsKeyPressed(KEY_BACKSPACE) && equationLetterCount > 0) {
                equationLetterCount--;
                equationText[equationLetterCount] = '\0'; // Clear the character at the current position
            } 

             if (key >= 32 && key <= 125) {
                equationText[equationLetterCount] = (char)key;
                equationLetterCount++;
            }
        }


        if (variableEditMode && equationLetterCount < 64) {
            int key = GetCharPressed();
                // Check for continuous backspace key press
            if (IsKeyDown(KEY_BACKSPACE)) {
                backspaceTimer += GetFrameTime();
                if (!backspaceHeld && backspaceTimer >= BACKSPACE_REPEAT_DELAY) {
                    backspaceHeld = true;
                    backspaceTimer = 0.0f;

                    if (varibaleLetterCount > 0) {
                        varibaleLetterCount--;
                        variableText[varibaleLetterCount] = '\0'; // Clear the character at the current position
                    }
                }
                else if (backspaceHeld && backspaceTimer >= BACKSPACE_REPEAT_INTERVAL) {
                    backspaceTimer = 0.0f;

                    if (varibaleLetterCount > 0) {
                        varibaleLetterCount--;
                        variableText[varibaleLetterCount] = '\0'; // Clear the character at the current position
                    }
                }
            }
            else {
                backspaceHeld = false;
                backspaceTimer = 0.0f;
            }

            if (IsKeyPressed(KEY_BACKSPACE) && varibaleLetterCount > 0) {
                varibaleLetterCount--;
                variableText[varibaleLetterCount] = '\0'; // Clear the character at the current position
            } 

            
            if ((key >= 48 && key <= 57) || key == 46) { // Only allow numbers and .
                variableText[varibaleLetterCount] = (char)key;
                varibaleLetterCount++;
                variableFloat = std::stof(variableText);
            }
            
        }

        BeginDrawing();
        ClearBackground(RAYWHITE);
        
        DrawText(title.c_str(), xPosition, yPosition/5, titleSize, BLACK);
        
        DrawText("Input Field:", equtionInputBox.x, equtionInputBox.y-equtionInputBoxHeight, 20, DARKGRAY);
        DrawRectangleLinesEx(equtionInputBox, 2, DARKGRAY);

        if (equationEditMode) {
            DrawRectangle(equtionInputBox.x + 1, equtionInputBox.y + 1, equtionInputBox.width - 2, equtionInputBox.height - 2, LIGHTGRAY);
        }

        DrawText(equationText, equtionInputBox.x + 10, equtionInputBox.y + 10, 20, MAROON);

        DrawText("Variable Field:", variableInputBox.x, variableInputBox.y-variableInputBoxHeight, 20, DARKGRAY);
        DrawRectangleLinesEx(variableInputBox, 2, DARKGRAY);

        if (variableEditMode) {
            DrawRectangle(variableInputBox.x + 1, variableInputBox.y + 1, variableInputBox.width - 2, variableInputBox.height - 2, LIGHTGRAY);
        }
        DrawText(variableText, variableInputBox.x + 10, variableInputBox.y + 10, 20, MAROON);
        // Calculate button
        const int calculateButtonWidth = 200;
        const int calculateButtonHeight = 40;
        Rectangle calculateButton = { (screen_width - calculateButtonWidth) / 2, (screen_height - calculateButtonHeight) / 1.5, calculateButtonWidth, calculateButtonHeight };
        bool calculateButtonClicked = false;

        if (CheckCollisionPointRec(GetMousePosition(), calculateButton) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            calculateButtonClicked = true;
        }

        if (calculateButtonClicked) {
            get_derivative(equationText, variableFloat);
            
        }

        // Draw calculate button
        DrawRectangleLinesEx(calculateButton, 2, DARKGRAY);
        DrawText("Calculate", calculateButton.x + 10, calculateButton.y + 10, 20, MAROON);

        //update the screen
        // Read derivative values from output.txt
        std::ifstream inputFile("output.txt");
        if (inputFile.is_open()) {
            std::string derivativeValue;
            int derivativeCount = 1;
            std::getline(inputFile, derivativeValue); // Skip the first line
            while (std::getline(inputFile, derivativeValue) && derivativeValue != "0") {
            // Show derivative value
            DrawText((std::to_string(derivativeCount) + "  derivative:  " + derivativeValue).c_str(), xPosition-100, calculateButton.y + 10 + (derivativeCount * 20), 20, MAROON);
            derivativeCount++;
            }

            inputFile.close();
        }

        inputFile.close();

        EndDrawing();
    }
    CloseWindow();
    return 0;
}
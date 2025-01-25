#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define RPN_CAPACITY 10

typedef struct stack_t {
    int top; //Index of the top stack value
    int data[RPN_CAPACITY]; //Stack data
} stack_t;

int is_full(stack_t* stack) {
    if(stack->top == RPN_CAPACITY) {
        return 1;
    } else {
        return 0;
    }
}

int is_empty(stack_t* stack) {
    if(stack->top == 0) {
        return 1;
    } else {
        return 0;
    }
}

void stack_init(stack_t *stack) {
    stack->top = 0;
}

void push(stack_t* stack, int num) {
    stack->top++;
    stack->data[stack->top] = num;
}

int pop(stack_t* stack) {
    int ret_num = stack->data[stack->top];
    stack->top--;
    return ret_num;
}

int main() {
    stack_t stack;
    int c, operand1, operand2, result;
    int line = 1;
    int overflow = 0;
    int not_enough_nums = 0;
    stack_init(&stack);

    while(1) {
        c = getchar();
        if(c == EOF) {
            break;
        } else if(overflow) {
            if(c == '\n') {
                overflow = 0;
                line++;
            }
            continue;
        } else if(not_enough_nums) {
            not_enough_nums = 0;
            line++;
            continue;
        } else if(isdigit(c)) {
            ungetc(c, stdin);
            int number;
            if(scanf("%d", &number) != 1) {
                printf("line %d: error at %c\n", line, c);
                line++;
                while((c = getchar()) != '\n' &&  c != EOF);
                continue;
            } else if(stack.top == RPN_CAPACITY) {
                printf("line %d: error at %d\n", line, number);
                while(!is_empty(&stack)) {
                    pop(&stack);
                }
                overflow = 1;
                continue;
            }
            push(&stack, number);
        } else if(c == '+' || c == '-' || c == '*' || c == '/') {
            if(stack.top < 2) {
                printf("line %d: error at %c\n", line, c);
                pop(&stack);
                not_enough_nums = 1;
                continue;
            } else {
                operand1 = pop(&stack);
                operand2 = pop(&stack);
                switch(c) {
                    case '+':
                        result = operand2 + operand1;
                        push(&stack, result);
                        break;
                    case '-':
                        result = operand2 - operand1;
                        push(&stack, result);
                        break;
                    case '*':
                        result = operand2 * operand1;
                        push(&stack, result);
                        break;
                    case '/':
                        if(operand1 == 0) {
                            printf("line %d: error at /\n", line);
                            line++;
                            while(!is_empty(&stack)) {
                                pop(&stack);
                            }
                            while((c = getchar()) != '\n' && c != EOF);
                            continue;
                        } else {
                            result = operand2 / operand1;
                            push(&stack, result);
                        }
                        break;
                    default:
                        break;

                }
            }
        } else if(c == '\n') {
            if(stack.top == 1) {
                printf("line %d: %d\n", line, pop(&stack));
            } else if(stack.top == 0) {
                printf("line %d: error at \\n\n", line);
            } else {
                printf("line %d: error at \\n\n", line);
                while(!is_empty(&stack)) {
                    pop(&stack);
                }
            }
            line++;
        } else if (isspace(c)) {
            continue;
        } else {
            printf("line %d: error at %c\n", line, c);
            while ((c = getchar()) != '\n' && c != EOF);
        }
    }
    return 0;
}


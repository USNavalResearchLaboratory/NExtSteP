#ifndef LINKEDlIST
#define LINKEDlIST 1

typedef struct llNodeStruct {
    void* contents;
    struct llNodeStruct* next;
} llNode;

typedef struct llRootStruct {
    int length;
    llNode* head;
    llNode* tail;
} llRoot;

/*
    Creates the root of a linked list.  The root contains a count of
    length of the list and pointers to the front and back of the
    list.

    Returns:
        llRoot*: the ponter to the created root.
*/
llRoot* createLinkedList(void);

/*
    Create a node of the linked list. A node contains a void* to some
    data and a pointer to the next item in the list.

    Input Args:
        void* value: The pointer you want inserted into the list.
    Returns:
        llNode*: The pointer to the created node.
*/
llNode* createllNode(void* value);

/*
    Checks to see if the linked list is empty.

    Input Aargs:
        llRoot* r: A pointer to the root of the list to check.
    Returns:
        int: 1 if true, 0 if false.
*/
int isEmpty(llRoot* r);

/*
    Insert node at the back of a list.

    Input Args:
        llRoot* r: The root of the linked list.
        void* value: The value to be inserted into the list.
*/
void insertllBack(llRoot* r, void* value);

/*
    Insert node at the front of a list.

    Input Args:
        llRoot* r: The root of the linked list.
        void* value: The value to be inserted into the list.
*/
void insertllFront(llRoot* r, void* value);

/*
    Frees the linked list. Only frees the linked list, anything
    pointed to by the list is not freed.

    Input Args:
        llRoot* r: the root of the linked list.
*/
void freeLL(llRoot* r);

/*
    Walks the list and passes the pointer to llNode->contents to a
    user supplied void function which takes pointer to a void* 
    (void**) an argument.

    Input Args:
        llRoot* r: The root of the linked list.
        void (*func)(): The 
*/
void applyToList(llRoot* r, void (*func)(void**));

/*
    A generic print function for the linked list. Uses a user defined
    function to print each item of the list.  The user defined
    function takes a pointer to the contents of the node (void**) of
    the linked list. If the pointer to the user defined funciton is
    NULL the pointer to each element of the list is printed instead.

    A newline is printed at the end of the list.

    If the list is empty the word "empty" followed by a newline char
    will be printed.

    Input Args:
        llRoot* r: The root of the linked list.
        void (*printFunction)(): A pointer to a user defined print
            function.

*/
void printList(llRoot* r, void (*printFunction)(void**));
#endif
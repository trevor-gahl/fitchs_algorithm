import sys


class TreeNode:
    def __init__(self, value, weight):
        self.left = None
        self.right = None
        self.data = value
        self.weight = weight


class Tree:
    def __init__(self):
        self.root = None

    def addNode(self, node, value, weight):
        if(node == None):
            self.root = TreeNode(value, weight)
        else:
            if(weight < node.weight):
                if(node.left == None):
                    node.left = TreeNode(value, weight)
                else:
                    self.addNode(node.left, value, weight)
            else:
                if(node.right == None):
                    node.right = TreeNode(value, weight)
                else:
                    self.addNode(node.right, value, weight)

    def printInorder(self, node):
        if(node != None):
            self.printInorder(node.left)
            print(node.data)
            self.printInorder(node.right)


def main():
    testTree = Tree()
    testTree.addNode(testTree.root, 'abc', 15)
    testTree.addNode(testTree.root, 'cba', 10)
    testTree.addNode(testTree.root, 'cca', 3)
    testTree.addNode(testTree.root, 'cbc', 8)
    testTree.printInorder(testTree.root)


if __name__ == '__main__':
    sys.exit(main())

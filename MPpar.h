#ifndef _PARSE_H
#define _PARSE_H

/* Port the parse.perl program */

#include <stdio.h>
#include "MPtag.h"

#define MAXBUFF 100000		// maximum length of input/output, ok to set very large

class MPparnode;

/*! \brief A class to perform parsing using simple part-of-speech regular expressions.
 *
 * This class inherits the MPtag class, and parsing automatically calls it to obtain
 * parts-of-speech. The regular expressions used to recognize phrases have not been
 * published or reviewed, and so there is little documentation for how to use or
 * interpret these results. See example code.
 */

class MPpar : public MPtag
{
public:
	friend class MPparnode;

	MPpar();				///< \brief A parse tree, can be reused to parse strings.
	MPpar(int nn, MPparnode **nd);		// Create a tree with the given nodes
	~MPpar();
	void copy_options(MPpar *t);		// Copy the options

	void load_str(const string& str);	// tag a string and make it the current tree, ready for parsing

	// Functions for parsing
	void parse(const string&);		///< \brief Load string and parse it
	void parse();				// Apply all parse rules to the current tree
	void match(const string& str, const string& pat, vector<string>& lbuff);	// Load string and do matching
	void match(const string& pat, vector<string>& lbuff);	// Extract all matching patterns from current tree
	void match(const string& str, const string& pat, vector<string>& lbuff, vector<string>& tags);	// Load string and do matching
	void match(const string& pat, vector<string>& lbuff, vector<string>& tags);	// Extract all matching patterns with tags from current tree
	void match(const string& pat, const string& tag);	// Replace all matching patterns with tag from current tree
	void parse_parens();			// Reorg and parse contents of parentheticals

	// Functions for printing

	void print(int rec);			///< \brief Print tree (with recurson)
	void copy_tokens(char *buff);		// Copy all of the tokens of the tree
	void copy_tags(char *buff);		// ... the tags
	void print_tokens();			// Print them out
	void print_nodes(const char *tg);
	void list_nodes(const char *tg, vector<string>& lbuff); ///< \brief After a parse, obtain a vector of constituents with tag \p tg
	void delete_nodes();			// Delete the nodes

	void flatten();				///< \brief Concatenate all nodes at top level of this tree
	void flatten(const char *tg);		///< \brief Concatenate all nodes anywhere in this tree having the given tag

	void untag();				// Convert parse tokens and tags to POS tags!?

	// These functions are used to perform attachment
	// of postmodifying phrases

	void transport(const char *tg, MPpar *t1, int n1, int n2);
	void search_right(int n, const char *tg, MPpar *&t1, int& n1);
	void attach(const char *pre, const char *post);
	void postmod();				/// \brief Perform all built-in attachments, e.g. NP to PP

	// Maintain ownership

	void set_owners();			// Set owner of all subtrees (needed for internal use only)

	char	*tagstr;			// The original text
	int	num_nodes;			// The number of nodes in the "current" parse
	MPparnode **nodes;			// An array of (pointers to) the nodes of the current parse
	MPpar	*owner;				// The owner of this tree (if any)

	int option_parens;			///< \brief Parse parens (default 0, 1=as is, 2=as constituent)
	int option_conj;			///< \brief Process coords (default 0, 1=simple, 2=phrase)
	int option_obj;				///< \brief Attach objects and arguments of verb phrases (defualt 0)
	int option_postmod;			///< \brief Attach postmodifying phrases (PP, VPP, etc) (default 0)
	int option_subparse;			///< \brief Parse deleted (parenthetical) subtrees (default 0)
	int option_printnull;			///< \brief Print removed trees
	int option_prep;			///< \brief Match prepositional phrases (default 1)

	int option_rev;				// Work right to left in string (internal only, default 0)
};

#define NODE_NONE 0
#define NODE_TOKEN 1
#define NODE_TREE 2

class MPparnode
{
public:
	MPparnode(const char *tg, const char *tk);	// Make a TOKEN node
	MPparnode(const char *tg, MPpar *tr);	// Make a TREE node
	~MPparnode();

	void flatten();				// Flatten a tree node into a token node

	int	type;
	char	*tag;
	char	*tok;
	MPpar	*tree;
};

#endif

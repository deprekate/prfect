--
-- Recode2 public schema 1.2
--  http://recode.ucc.ie
--  copyright 2009
--

SET client_encoding = 'UTF8';
SET check_function_bodies = false;
SET client_min_messages = warning;

SET search_path = public, pg_catalog;

--
-- Name: SCHEMA public; Type: COMMENT; Schema: -; Owner: postgres
--

COMMENT ON SCHEMA public IS 'Recode2 Public Database schema';


SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: events; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE events (
    event integer NOT NULL,
    name character varying(32) NOT NULL,
    description text
);


--
-- Name: journals; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE journals (
    id serial NOT NULL,
    id_recode character varying(64) NOT NULL,
    authors text NOT NULL,
    title text NOT NULL,
    journal character varying(128) NOT NULL,
    theyear integer NOT NULL,
    doi character varying(64)
);


--
-- Name: kingdoms; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE kingdoms (
    kingdom integer NOT NULL,
    phylum character varying(32) NOT NULL,
    genome character varying(32) NOT NULL,
    description text
);


--
-- Name: molecules; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE molecules (
    id serial NOT NULL,
    id_recode character varying(64) NOT NULL,
    variant character varying(128),
    ext_databases text,
    coordinates character varying(64) DEFAULT NULL::character varying,
    "sequence" text NOT NULL,
    genetic integer DEFAULT 1 NOT NULL,
    acid_type character varying(8),
    strand integer DEFAULT 0 NOT NULL,
    length bigint DEFAULT 0 NOT NULL,
    topology character varying(8),
    comments text,
    annotation text
);


--
-- Name: morphorna; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE morphorna (
    id serial NOT NULL,
    id_recode character varying(64) NOT NULL,
    id_molecule integer NOT NULL,
    structure integer DEFAULT 0 NOT NULL,
    description text NOT NULL
);


--
-- Name: organisms; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE organisms (
    latin character varying(128) NOT NULL,
    variant character varying(128),
    acronym character varying(32),
    synonym text,
    taxonid bigint NOT NULL,
    genus character varying(128),
    phylum text
);


--
-- Name: products; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE products (
    id serial NOT NULL,
    id_product character varying(64) NOT NULL,
    id_molecule integer NOT NULL,
    id_recode character varying(64) NOT NULL,
    name character varying(64) NOT NULL,
    "sequence" text,
    modification text,
    coordinates text,
    description text
);


--
-- Name: recode2; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE recode2 (
    id serial NOT NULL,
    id_recode character varying(64) NOT NULL,
    legacy character varying(128) DEFAULT NULL::character varying,
    organism bigint NOT NULL,
    kingdom integer DEFAULT 0 NOT NULL,
    locus character varying(64) NOT NULL,
    description text,
    status integer DEFAULT 0
);


--
-- Name: recoding; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE recoding (
    id serial NOT NULL,
    id_recode character varying(64) NOT NULL,
    id_product character varying(64) NOT NULL,
    event integer DEFAULT 0 NOT NULL,
    experimental character varying(1) DEFAULT 'f'::character varying NOT NULL,
    "position" bigint,
    codon character varying(64),
    upstream character varying(64),
    esite character varying(64),
    asite character varying(64),
    psite character varying(64),
    downstream character varying(64),
    model text,
    description text
);


--
-- Name: status; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE status (
    status integer NOT NULL,
    name character varying(32) NOT NULL,
    description text
);


--
-- Name: structures; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE structures (
    structure integer NOT NULL,
    name character varying(32) NOT NULL,
    description text
);


--
-- Name: translation; Type: TABLE; Schema: public; Owner: apache; Tablespace: 
--

CREATE TABLE translation (
    id integer NOT NULL,
    translation character varying(128) NOT NULL
);


--
-- Data for Name: events; Type: TABLE DATA; Schema: public; Owner: apache
--

INSERT INTO events (event, name, description) VALUES (0, 'Unknown', 'The recoding event is unknown or not classified');
INSERT INTO events (event, name, description) VALUES (1, '-1 Frameshifting', 'A directed change in <a href="http://en.wikipedia.org/wiki/Translation_%28genetics%29">translational</a> reading frame in -1 that allows the production of two proteins from a single <a href="http://en.wikipedia.org/wiki/Messenger_RNA">mRNA</a>. The process is programmed by the nucleotide sequence of the mRNA and is sometimes also affected by the secondary or tertiary mRNA <a href="http://en.wikipedia.org/wiki/Secondary_structure">structure</a>');
INSERT INTO events (event, name, description) VALUES (2, '+1 Frameshifting', 'A directed change in <a href="http://en.wikipedia.org/wiki/Translation_%28genetics%29">translational</a> reading frame in +1 that allows the production of two proteins from a single <a href="http://en.wikipedia.org/wiki/Messenger_RNA">mRNA</a>. The process is programmed by the nucleotide sequence of the mRNA');
INSERT INTO events (event, name, description) VALUES (3, 'Readthrough', 'A <a href="http://en.wikipedia.org/wiki/Codon">codon</a> redefinition occurs when a stop codon specifies a standard <a href="http://en.wikipedia.org/wiki/Amino_acid">amino acid</a> (often <a href="http://en.wikipedia.org/wiki/Glutamine">glutamine</a> or <a href="http://en.wikipedia.org/wiki/Tryptophan">tryptophan</a>)');
INSERT INTO events (event, name, description) VALUES (4, 'Selenocysteine', 'A <a href="http://en.wikipedia.org/wiki/Codon">codon</a> redefinition occurs when a stop codon specifies a <a href="http://en.wikipedia.org/wiki/Selenocysteine">selenocysteine</a> insertion');
INSERT INTO events (event, name, description) VALUES (5, 'Pyrolysine', 'A <a href="http://en.wikipedia.org/wiki/Codon">codon</a> redefinition occurs when a stop codon specifies a pyrolysine insertion');
INSERT INTO events (event, name, description) VALUES (6, 'Hopping', 'A <a href="http://en.wikipedia.org/wiki/Translation_%28genetics%29">translational</a> bypass allows the coupling of two <a href="http://en.wikipedia.org/wiki/ORF">ORFs</a> separated on an <a href="http://en.wikipedia.org/wiki/Messenger_RNA">mRNA</a> by a coding gap');


--
-- Data for Name: kingdoms; Type: TABLE DATA; Schema: public; Owner: apache
--

INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (0, 'Unknown', 'Unknown', 'The recoding is not classified');
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (1, 'Archaea', 'Chromosome', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (2, 'Archaea', 'Plasmid', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (3, 'Archaea', 'Insertion sequence', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (4, 'Bacteria', 'Chromosome', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (5, 'Bacteria', 'Plasmid', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (6, 'Bacteria', 'Insertion sequence', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (7, 'Bacteria', 'Prophage', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (8, 'Eukaryota', 'Chromosome', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (9, 'Eukaryota', 'Plasmid', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (10, 'Eukaryota', 'Transposable element', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (11, 'Eukaryota', 'Organelle', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (14, 'Viroid', 'Viroid', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (15, 'Synthetic', 'Synthetic', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (12, 'Viral Genome', 'Virus', NULL);
INSERT INTO kingdoms (kingdom, phylum, genome, description) VALUES (13, 'Viral Genome', 'Phage', NULL);


--
-- Data for Name: status; Type: TABLE DATA; Schema: public; Owner: apache
--

INSERT INTO status (status, name, description) VALUES (0, 'Unevaluated', 'The Recode record is not evaluated');
INSERT INTO status (status, name, description) VALUES (1, 'Inferred', 'The Recode record is automatically <a href="http://en.wikipedia.org/wiki/Inference">inferred</a> by <a href="http://en.wikipedia.org/wiki/Genome">genome</a> sequence analysis. There is no experimental support for the full extent of the recoding; there may be some level of support by <a href="http://en.wikipedia.org/wiki/Homology_%28biology%29">homology</a>');
INSERT INTO status (status, name, description) VALUES (2, 'Model', 'The Recode record is automatically <a href="http://en.wikipedia.org/wiki/Predicted">predicted</a> by <a href="http://en.wikipedia.org/wiki/Genome">genome</a> sequence analysis using a <a href="http://en.wikipedia.org/wiki/Model">Model</a> of recoding event. The record may represent an <em>ab initio</em> prediction, and may have some or <a href="http://en.wikipedia.org/wiki/Homology_%28biology%29">homology</a> support');
INSERT INTO status (status, name, description) VALUES (3, 'Predicted', 'The Recode record is automatically <a href="http://en.wikipedia.org/wiki/Predicted">predicted</a> and has not been subject to individual <a href="http://en.wikipedia.org/wiki/Review">review</a>. Strong similarity and <a href="http://en.wikipedia.org/wiki/Secondary_structure">structure</a> are supported');
INSERT INTO status (status, name, description) VALUES (4, 'Provisional', 'The Recode record has not yet been subject to <a href="http://en.wikipedia.org/wiki/Review">review</a>/<a href="http://en.wikipedia.org/wiki/Scientific_paper">publication</a>, but has been manually submit by outside collaborators');
INSERT INTO status (status, name, description) VALUES (5, 'Reviewed', 'The Recode record has been the <a href="http://en.wikipedia.org/wiki/Review">reviewed</a> by staff or is <a href="http://en.wikipedia.org/wiki/Scientific_paper">published</a>. The review process includes reviewing available sequence data and frequently also includes a review of the literature and other sources of information. Some Recode records may incorporate expanded sequence and annotation information including additional publications and features, as deemed relevant');
INSERT INTO status (status, name, description) VALUES (6, 'Validated', 'The Recode record is experimentally validated. The record is used as reference for <a href="http://en.wikipedia.org/wiki/Model">Model</a>, and similarity <a href="http://en.wikipedia.org/wiki/Inference">inferring</a>');
INSERT INTO status (status, name, description) VALUES (10, 'Removed', 'The Recode record was removed from the database. The record is still in the database itselft, but is not publicly avalable');


--
-- Data for Name: structures; Type: TABLE DATA; Schema: public; Owner: apache
--

INSERT INTO structures (structure, name, description) VALUES (0, 'Unknown', '');
INSERT INTO structures (structure, name, description) VALUES (1, 'Trans/Interaction', '');
INSERT INTO structures (structure, name, description) VALUES (2, 'Sequence', '');
INSERT INTO structures (structure, name, description) VALUES (3, 'Stemloop', '');
INSERT INTO structures (structure, name, description) VALUES (4, 'SECIS', '');
INSERT INTO structures (structure, name, description) VALUES (5, 'Pseudoknot', '');
INSERT INTO structures (structure, name, description) VALUES (6, 'Complex structure', 'Vienna format');
INSERT INTO structures (structure, name, description) VALUES (7, 'SECIS structure', 'Vienna (pseudoknot enable) format');
INSERT INTO structures (structure, name, description) VALUES (8, 'Pseudoknot structure', 'Vienna (pseudoknot enable) format');


--
-- Data for Name: translation; Type: TABLE DATA; Schema: public; Owner: apache
--

INSERT INTO translation (id, translation) VALUES (0, 'Unknown');
INSERT INTO translation (id, translation) VALUES (1, 'Standard');
INSERT INTO translation (id, translation) VALUES (2, 'Vertebrate Mitochondrial');
INSERT INTO translation (id, translation) VALUES (3, 'Yeast Mitochondrial');
INSERT INTO translation (id, translation) VALUES (4, 'Mold, Protozoan, Coelenterate Mito. and Myco/Spiroplasma');
INSERT INTO translation (id, translation) VALUES (5, 'Invertebrate Mitochondrial');
INSERT INTO translation (id, translation) VALUES (6, 'Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear');
INSERT INTO translation (id, translation) VALUES (9, 'Echinoderm Mitochondrial');
INSERT INTO translation (id, translation) VALUES (10, 'Euploid Nuclear');
INSERT INTO translation (id, translation) VALUES (11, 'Bacterial');
INSERT INTO translation (id, translation) VALUES (12, 'Alternative Yeast Nuclear');
INSERT INTO translation (id, translation) VALUES (13, 'Ascidian Mitochondrial');
INSERT INTO translation (id, translation) VALUES (14, 'Flatworm Mitochondrial');
INSERT INTO translation (id, translation) VALUES (15, 'Blepharisma Macronuclear');
INSERT INTO translation (id, translation) VALUES (16, 'Chlorophycean Mitochondrial');
INSERT INTO translation (id, translation) VALUES (21, 'Trematode Mitochondrial');


--
-- Name: events_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY events
    ADD CONSTRAINT events_pkey PRIMARY KEY (event);


--
-- Name: events_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY events
    ADD CONSTRAINT events_unique UNIQUE (name);


--
-- Name: journals_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY journals
    ADD CONSTRAINT journals_pkey PRIMARY KEY (id);


--
-- Name: kingdoms_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY kingdoms
    ADD CONSTRAINT kingdoms_pkey PRIMARY KEY (kingdom);


--
-- Name: molecules_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY molecules
    ADD CONSTRAINT molecules_pkey PRIMARY KEY (id);


--
-- Name: morphorna_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY morphorna
    ADD CONSTRAINT morphorna_pkey PRIMARY KEY (id);


--
-- Name: organisms_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY organisms
    ADD CONSTRAINT organisms_pkey PRIMARY KEY (taxonid);


--
-- Name: organisms_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY organisms
    ADD CONSTRAINT organisms_unique UNIQUE (latin, variant);


--
-- Name: products_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY products
    ADD CONSTRAINT products_pkey PRIMARY KEY (id);


--
-- Name: products_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY products
    ADD CONSTRAINT products_unique UNIQUE (id_product);


--
-- Name: recode2_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY recode2
    ADD CONSTRAINT recode2_pkey PRIMARY KEY (id);


--
-- Name: recode2_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY recode2
    ADD CONSTRAINT recode2_unique UNIQUE (id_recode);


--
-- Name: recoding_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY recoding
    ADD CONSTRAINT recoding_pkey PRIMARY KEY (id);


--
-- Name: status_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY status
    ADD CONSTRAINT status_pkey PRIMARY KEY (status);


--
-- Name: status_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY status
    ADD CONSTRAINT status_unique UNIQUE (name);


--
-- Name: structures_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY structures
    ADD CONSTRAINT structures_pkey PRIMARY KEY (structure);


--
-- Name: structures_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY structures
    ADD CONSTRAINT structures_unique UNIQUE (name);


--
-- Name: translation_pkey; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY translation
    ADD CONSTRAINT translation_pkey PRIMARY KEY (id);


--
-- Name: translation_unique; Type: CONSTRAINT; Schema: public; Owner: apache; Tablespace: 
--

ALTER TABLE ONLY translation
    ADD CONSTRAINT translation_unique UNIQUE (translation);


--
-- Name: journals_id_recode_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY journals
    ADD CONSTRAINT journals_id_recode_fkey FOREIGN KEY (id_recode) REFERENCES recode2(id_recode) ON DELETE CASCADE;


--
-- Name: molecules_genetic_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY molecules
    ADD CONSTRAINT molecules_genetic_fkey FOREIGN KEY (genetic) REFERENCES translation(id) ON DELETE CASCADE;


--
-- Name: molecules_id_recode_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY molecules
    ADD CONSTRAINT molecules_id_recode_fkey FOREIGN KEY (id_recode) REFERENCES recode2(id_recode) ON DELETE CASCADE;


--
-- Name: morphorna_id_recode_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY morphorna
    ADD CONSTRAINT morphorna_id_recode_fkey FOREIGN KEY (id_recode) REFERENCES recode2(id_recode) ON DELETE CASCADE;


--
-- Name: morphorna_structure_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY morphorna
    ADD CONSTRAINT morphorna_structure_fkey FOREIGN KEY (structure) REFERENCES structures(structure) ON DELETE CASCADE;


--
-- Name: products_id_recode_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY products
    ADD CONSTRAINT products_id_recode_fkey FOREIGN KEY (id_recode) REFERENCES recode2(id_recode) ON DELETE CASCADE;


--
-- Name: recode2_kingdom_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY recode2
    ADD CONSTRAINT recode2_kingdom_fkey FOREIGN KEY (kingdom) REFERENCES kingdoms(kingdom) ON DELETE CASCADE;


--
-- Name: recode2_organism_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY recode2
    ADD CONSTRAINT recode2_organism_fkey FOREIGN KEY (organism) REFERENCES organisms(taxonid) ON DELETE CASCADE;


--
-- Name: recoding_event_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY recoding
    ADD CONSTRAINT recoding_event_fkey FOREIGN KEY (event) REFERENCES events(event) ON DELETE CASCADE;


--
-- Name: recoding_id_product_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY recoding
    ADD CONSTRAINT recoding_id_product_fkey FOREIGN KEY (id_product) REFERENCES products(id_product) ON DELETE CASCADE;


--
-- Name: recoding_id_recode_fkey; Type: FK CONSTRAINT; Schema: public; Owner: apache
--

ALTER TABLE ONLY recoding
    ADD CONSTRAINT recoding_id_recode_fkey FOREIGN KEY (id_recode) REFERENCES recode2(id_recode) ON DELETE CASCADE;

--

